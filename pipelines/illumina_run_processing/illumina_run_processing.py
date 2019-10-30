#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
from __future__ import print_function, division, unicode_literals, absolute_import
import argparse
import logging
import os
import re
import sys
import itertools
import subprocess
import xml.etree.ElementTree as Xml
import math
import subprocess
import csv
import json
from collections import OrderedDict

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs

from bfx.readset import parse_illumina_raw_readset_files
from bfx import bvatools
from bfx import picard
from bfx import tools
from bfx import run_processing_tools
from bfx import fastqc
from bfx import bash_cmd as bash

from pipelines import common

log = logging.getLogger(__name__)

class RunInfoRead(object):
    """ 
    Model of a read from the Illumina sequencer.
    Those attributes can be found in the RunInfo.xml file.
    """

    def __init__(self, number, nb_cycles, is_index):
        self._number = number
        self._nb_cycles = nb_cycles
        self._is_index = is_index

    @property
    def number(self):
        return self._number

    @property
    def nb_cycles(self):
        return self._nb_cycles

    @property
    def is_index(self):
        return self._is_index

class IlluminaRunProcessing(common.MUGQICPipeline):
    """
    Illumina Run Processing Pipeline
    ================================

    The standard MUGQIC Illumina Run Processing pipeline uses the Illumina bcl2fastq
    software to convert and demultiplex the base call files to fastq files. The
    pipeline runs some QCs on the raw data, on the fastq and on the alignment.

    Sample Sheets
    -------------

    The pipeline uses two input sample sheets. The first one is the standard Casava
    sheet, a csv file having the following columns (please refer to the Illumina
    Casava user guide):

    - `SampleID`
    - `FCID`
    - `SampleRef`
    - `Index`
    - `Description`
    - `Control`
    - `Recipe`
    - `Operator`
    - `SampleProject`

    Example:

        FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
        H84WNADXX,1,sample1_MPS0001,,TAAGGCGA-AGAGTAGA,,N,,,nanuq
        H84WNADXX,1,sample47_MPS0047,,GTAGAGGA-CTAAGCCT,,N,,,nanuq

    The second sample sheet is called the Nanuq run sheet. It's a csv file with the
    following minimal set of mandatory columns (the column order in the file doesn't
    matter)

    - `ProcessingSheetId` Must be the same as the `SampleID` from the Casava Sheet.
    - `Name` The sample name put in RG headers of bam files and on filename on disk.
    - `Run` The run number.
    - `Region` The lane number.
    - `Library Barcode` The library barcode put in .bam's RG headers and on disk
    - `Library Source` The type of library. If this value contains `RNA` or `cDNA`,
    `STAR` will be used to make the aligmnent, otherwise, `bwa_mem` will be used
    - `Library Type` Used to determine is the sample is from cDNA/RNA when the
    `Library Source` is `Library`
    - `BED Files` The name of the BED file containing the genomic targets. This is
    the `filename` parameter passed to the `fetch_bed_file_command`
    - `Genomic Database` The reference used to make the alignment and calculate aligments metrics

    Example:

        Name,Genomic Database,Library Barcode,Library Source,Library Type,Run,Region,BED Files,ProcessingSheetId
        sample1,Rattus_norvegicus:Rnor_5.0,MPS0001,RNA,Nextera XT,1419,1,toto.bed,sample1_MPS0001
        sample47,,MPS1047,Library,Nextera XT,1419,2,toto.bed,sample47_MPS1047
    """

    def __init__(self, protocol=None):
        self._protocol=protocol
        self.copy_job_inputs = {}
        self.argparser.add_argument("-d", "--run", help="Run directory (mandatory)", required=False, dest="run_dir")
        self.argparser.add_argument("--run-id", help="Run ID. Default is parsed from the run folder", required=False, dest="run_id")
#        self.argparser.add_argument("--flowcell-id", help="Flowcell ID. Default is parsed from the run folder", required=False, dest="flowcell_id")
        self.argparser.add_argument("--raw-fastq-prefix", help="Prefix used to search for the raw fastq from the sequencer. Default <FLOWCELL_ID>_<RUN_ID>", required=False, dest="raw_fastq_prefix")
        self.argparser.add_argument("--lane", help="Lane number (to only process the given lane)", type=int, required=False, dest="lane_number")
        self.argparser.add_argument("-r", "--readsets", help="Readset file i.e. Clarity event file (mandatory)", type=file, required=False)
        self.argparser.add_argument("-x", help="First index base to use for demultiplexing (inclusive). The index from the sample sheet will be adjusted according to that value.", type=int, required=False, dest="first_index")
        self.argparser.add_argument("-y", help="Last index base to use for demultiplexing (inclusive)", type=int, required=False, dest="last_index")
        self.argparser.add_argument("-m", help="Number of index mistmaches allowed for demultiplexing (default 1). Barcode collisions are always checked.", type=int, required=False, dest="number_of_mismatches")

        super(IlluminaRunProcessing, self).__init__(protocol)

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = {}
            if not hasattr(self, "_mask"):
                self._mask = {}
            for lane in self.lanes:
                self._readsets[lane] = self.load_readsets(lane)
                self._mask[lane] = self.get_mask(lane)
                self.generate_illumina_lane_casava_sheet(lane)
        return self._readsets

    @property
    def samples(self):
        if not hasattr(self, "_samples"):
            self._samples = {}
        for lane in self.lanes:
            self._samples[lane] = list(OrderedDict.fromkeys([readset.sample for readset in self.readsets[lane]]))
        return self._samples

    @property
    def lanes(self):
        if not hasattr(self, "_lanes"):
            if self.lane_number:
                self._lanes = [lane]
            else:
                self._lanes = [lane for lane in list(set([line['Position'].split(":")[0] for line in csv.DictReader(open(self.readset_file, 'rb'), delimiter=str(u'\t').encode('utf-8'), quotechar=str(u'"').encode('utf-8'))]))]
            for lane in self._lanes:
                self.copy_job_inputs[lane] = []
        return self._lanes

    @property
    def is_paired_end(self):
        if not hasattr(self, "_is_paired_end"):
            if len([read_info for read_info in self.read_infos if not read_info.is_index]) > 1:
                self._is_paired_end = True
            else:
                self._is_paired_end = False
        return self._is_paired_end

    @property
    def is_dual_index(self):
        if not hasattr(self, "_is_dual_index"):
            if len([read_info for read_info in self.read_infos[lane] if read_info.is_index]) > 1:
                self._is_dual_index = True
            else:
                self._is_dual_index = False
        return self._is_dual_index

    @property
    def mask(self):
        if not hasattr(self, "_mask"):
            _raise(SanitycheckError("No mask could be found !!"))
        return self._mask 

    @property
    def merge_undetermined(self):
        if not hasattr(self, "_merge_undetermined"):
            self._merge_undetermined = {}
            for lane in self.lanes:
                self._merge_undetermined[lane] = False
                # If only one library on the lane
                if len(self.readsets[lane]) == 1:
                    self._merge_undetermined[lane] = True
                if config.param('fastq', 'merge_undetermined', required=False, type='boolean'):
                    self._merge_undetermined[lane] = config.param('fastq', 'merge_undetermined')
        return self._merge_undetermined

    @property
    def instrument(self):
        if not hasattr(self, "_instrument"):
            self._instrument = ""
            for lane in self.lanes:
                lane_instrument = self.get_instrument()
                if self._instrument and self._instrument != lane_instrument:
                    _raise(SanitycheckError("One run (\"" + self.run_id + "\") cannot be shared by different instruments !! (\"" + lane_instrument + "\" vs. \"" +  self._instrument + "\")"))
                else:
                    self._instrument = lane_instrument
        return self._instrument

    @property
    def experiment_name(self):
        if not hasattr(self, "_experiment_name"):
            self._experiment_name = self.get_experiment_name()
        return self._experiment_name

    @property
    def run_number(self):
        if not hasattr(self, "_run_number"):
            self._run_number = self.get_run_number()
        return self._run_number

    @property
    def seq_category(self):
        if not hasattr(self, "_seq_category"):
            self._seq_category = self.get_seq_category()
        return self._seq_category

    @property
    def run_id(self):
        """
        The run id from the run folder.
        Supports both default folder name configuration and GQ's globaly unique name convention.
        """
        if not hasattr(self, "_run_id"):
            if re.search(".*_\d+HS\d\d[AB]", self.run_dir):
                m = re.search(".*/(\d+_[^_]+_\d+_[^_]+_(\d+)HS.+)", self.run_dir)
                self._run_id = m.group(2)
            elif re.search(".*\d+_[^_]+_\d+_.+", self.run_dir):
                m = re.search(".*/(\d+_([^_]+_\d+)_.*)", self.run_dir)
                self._run_id = m.group(2)
            else:
                log.warn("Unsupported folder name: " + self.run_dir)

        return self._run_id

    @property
    def run_dir(self):
        if self.args.run_dir:
            tmp_dir = os.path.expandvars(self.args.run_dir)
            if not os.path.isabs(tmp_dir):
                # Get relative path from output directory
                tmp_dir = os.path.relpath(tmp_dir, self.output_dir)
                # And get the absolute path
                tmp_dir = os.path.normpath(os.path.join(self.output_dir, tmp_dir))
            return tmp_dir
        else:
            _raise(SanitycheckError("Error: missing '-d/--run' option!"))

    @property
    def flowcell_id(self):
        if not hasattr(self, "_flowcell_id"):
            self._flowcell_id = ""
            for lane in self.lanes:
                fc_set = set([readset.flow_cell for readset in self.readsets[lane] if readset.flow_cell])
                if len(fc_set) == 0:
                    _raise(SanitycheckError("Error: no flowcell id could be found in the readset objects for lane '" + lane + "'..."))
                if len(fc_set) > 1:
                    _raise(SanitycheckError("Error: more than one flowcell id found in the readset objects for lane '" + lane + "'..."))
                fc_id = list(fc_set)[0]
                if not self._flowcell_id:
                    self._flowcell_id = fc_id
                elif not self._flowcell_id == fc_id:
                    _raise(SanitycheckError("Error: more than one flowcell id found for the run, please check the readset file '" + self.readset_file + "'...")) 
        return self._flowcell_id

    @property
    def lane_number(self):
        if self.args.lane_number:
            return self.args.lane_number
        else:
            None

    @property
    def readset_file(self):
        if not hasattr(self, "_readset_file"):
            if self.args.readsets:
                self._readset_file = os.path.realpath(self.args.readsets.name)
            else:
                _raise(SanitycheckError("Error: missing '-r/--readsets' argument !"))
        return self._readset_file

    @property
    def bcl2fastq_job_input(self):
        if self.protocol == "clarity":
            return self.clarity_event_file
        else:
            return self.casava_sheet_file

    @property
    def number_of_mismatches(self):
        return self.args.number_of_mismatches if (self.args.number_of_mismatches is not None) else 1

    @property
    def first_index(self):
        return self.args.first_index if self.args.first_index else 1

    @property
    def last_index(self):
        return self.args.last_index if self.args.last_index else 999

    @property
    def umi(self):
        if not hasattr(self, "_umi"):
            self._umi = False
        return self._umi

    @property
    def bcl2fastq_extra_option(self):
        if not hasattr(self, "_bcl2fastq_extra_option"):
            return ""
        return self._bcl2fastq_extra_option

    @property
    def index1cycles(self):
        if not hasattr(self, "_index1cycles"):
            [self._index1cycles, self._index2cycles] = self.get_indexcycles()
        return self._index1cycles

    @property
    def index2cycles(self):
        if not hasattr(self, "_index2cycles"):
            [self._index1cycles, self._index2cycles] = self.get_indexcycles()
        return self._index2cycles

    @property
    def index_per_readset(self):
        if not hasattr(self, "_index_per_readset"):
            return ""
        # Define in generate_clarity_sample_sheet() 
        return self._index_per_readset

    @property
    def seqtype(self):
        if not hasattr(self, "_seqtype"):
            self._seqtype = self.get_seqtype()
        return self._seqtype

    @property
    def read_infos(self):
        if not hasattr(self, "_read_infos"):
            self._read_infos = self.parse_run_info_file()
        return self._read_infos

    @property
    def report_hash(self):
        if not hasattr(self, "_report_hash"):
            self._report_hash = {
                "run" : self.run_id,
                "instrument" : self.instrument,
                "flowcell" : self.flowcell_id,
                "lane" : self.lane_number,
                "seqtype" : self.seqtype,
                "sequencing_method" : "PAIRED_END" if self.is_paired_end else "SINGLE_END",
                "steps" : [],
                "run_validation" : [] 
            }
        return self._report_hash

    @property
    def report_inputs(self):
        if not hasattr(self, "_report_inputs"):
            self._report_inputs = {}
            for lane in self.lanes:
                self._report_inputs[lane] = {
                    'index' : {},
                    'sample_tag' : {},
                    'qc' : {},
                    'blast' : {},
                    'mark_dup' : {},
                    'align' : {}
                }
        return self._report_inputs

    def index(self):
        """
        Generate a file with all the indexes found in the index-reads of the run.

        The input barcode file is a two columns tsv file. Each line has a
        `barcode_sequence` and the corresponding `barcode_name`. This file can be
        generated by a LIMS.

        The output is a tsv file named `RUNFOLDER_LANENUMBER.metrics` that will be
        saved in the output directory. This file has four columns, the barcode/index
        sequence, the index name, the number of reads and the number of reads that have
        passed the filter.
        """
        jobs = []

        mask = ""
        index_length = self.get_sequencer_index_length()

        for read in self.read_infos:
            if read.is_index:
                mask += str(index_length) + 'B'
                break
            else:
                mask += str(read.nb_cycles) + 'T'

        if index_length == 0:
            log.info("No Indexes, *NOT* Generating index counts")

        else:
            for lane in self.lanes:
                lane_jobs = []

                input = os.path.join(self.run_dir, "RunInfo.xml")
                output = os.path.join(self.output_dir, "index",  self.run_id + "_" + lane + '.metrics')
                basecalls_dir = os.path.join(self.output_dir, "index", "BaseCalls")

                # CountIlluminaBarcode
                lane_jobs.append(
                    concat_jobs([
                        bash.mkdir(basecalls_dir),
                        bash.ln(
                            os.path.join(self.run_dir, "Data", "Intensities", "BaseCalls", "L00" + lane),
                            os.path.join(basecalls_dir, "L00" + lane)
                        ),
                        bash.ln(
                            os.path.join(self.run_dir, "Data", "Intensities", "s.locs"),
                            os.path.join(self.output_dir, "index", "s.locs")
                        ),
                        run_processing_tools.index(
                            input,
                            output,
                            basecalls_dir,
                            self.number_of_mismatches,
                            lane,
                            mask
                    )],
                    name="index." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                ))

                # BCL2fasq
                lane_jobs.append(
                    concat_jobs([
                        run_processing_tools.bcl2fastq_for_index(
                            self.run_dir,
                            os.path.join(self.output_dir, "index", "L00" + lane),
                            os.path.join(self.output_dir, "casavasheet." + lane + ".indexed.csv"),
                            self.flowcell_id,
                            lane,
                            demultiplex=True,
                            mismatches=self.number_of_mismatches,
                            mask=self.mask[lane]
                        ),
                        bash.cp(
                            os.path.join(self.output_dir, "index", "L00" + lane, "Reports", "html", self.flowcell_id, "all/all/all/lane.html"),
                            os.path.join(self.output_dir, "index", self.run_id + "_" + lane + ".index_stats.html")
                        ),
                        tools.edit_index_stats(
                            os.path.join(self.output_dir, "index", self.run_id + "_" + lane + ".index_per_readset.json"),
                            os.path.join(self.output_dir, "index", "L00" + lane, "Stats", "Stats.json"),
                            os.path.join(self.output_dir, "index", self.run_id + "_" + lane + ".index_stats.json")
                        )
                    ],
                    name="bcl2fastq_index." + self.run_id + "." + lane,
                    samples=self.samples[lane],
                    removable_files=[os.path.join(self.output_dir, "index", "L00" + lane)]
                ))

                # Index Validation
                idx_val_job = tools.index_validation(
                    os.path.join(self.output_dir, "index", self.run_id + "_" + lane + ".index_per_readset.json"),
                    output,
                    os.path.join(self.output_dir, "index", self.run_id + "_" + lane + ".index_stats.html"),
                    lane,
                    self.number_of_mismatches
                )
                idx_val_job.name = "index_validation." + self.run_id + "." + lane
                idx_val_job.samples = self.samples[lane]
                lane_jobs.append(idx_val_job)

                self.add_to_report_hash("index", lane_jobs)
                for readset in self.readsets[lane]:
                    self.report_inputs[lane]['index'][readset.name] = [
                        os.path.join(self.output_dir, "index", self.run_id + "_" + lane + ".index_stats.json"),
                        os.path.join(self.output_dir, "index", self.run_id + "_" + lane + ".index_per_readset.json")
                    ]
                self.add_copy_job_inputs(lane_jobs, lane)
                jobs.extend(lane_jobs)

        return jobs

    def fastq(self):
        """
        Launch fastq generation from Illumina raw data using BCL2FASTQ conversion
        software.

        The index base mask is calculated according to the sample and run configuration;
        and also according the mask parameters received (first/last index bases). The
        Casava sample sheet is generated with this mask. The default number of
        mismatches allowed in the index sequence is 1 and can be overrided with an
        command line argument. Demultiplexing always occurs even when there is only one
        sample in the lane, because then we merge undertermined reads.

        An optional notification command can be launched to notify the start of the
        fastq generation with the calculated mask.
        """
        jobs = []

        for lane in self.lanes:
            lane_jobs = []

            input = self.readset_file
            fastq_outputs, final_fastq_jobs = self.generate_fastq_outputs(lane)
            output_dir = os.path.join(self.output_dir, "Unaligned." + lane)
            casava_sample_sheet = os.path.join(self.output_dir, "casavasheet." + lane + ".indexed.csv")

            if self.umi:
                output_dir_noindex = os.path.join(self.output_dir, "Unaligned." + lane + ".noindex")
                casava_sample_sheet_noindex = os.path.join(self.output_dir, "casavasheet." + lane + ".noindex.csv")

                lane_jobs.append(
                    concat_jobs([
                        run_processing_tools.bcl2fastq(
                            input,
                            fastq_outputs,
                            output_dir,
                            casava_sample_sheet,
                            self.run_dir,
                            lane,
                            self.bcl2fastq_extra_option,
                            demultiplex=True,
                            mismatches=self.number_of_mismatches,
                            mask=self.mask[lane]
                        ),
                        run_processing_tools.bcl2fastq(
                            input,
                            fastq_outputs,
                            output_dir_noindex,
                            casava_sample_sheet_noindex,
                            self.run_dir,
                            lane,
                            self.bcl2fastq_extra_option
                    )],
                    name="fastq." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                ))

            else:
                bcl2fastq_job = run_processing_tools.bcl2fastq(
                    input,
                    fastq_outputs,
                    output_dir,
                    casava_sample_sheet,
                    self.run_dir,
                    self.lane_number,
                    self.bcl2fastq_extra_option,
                    demultiplex=True,
                    mismatches=self.number_of_mismatches,
                    mask=self.mask[lane]
                )
                bcl2fastq_job.name = "fastq." + self.run_id + "." + lane
                bcl2fastq_job.samples = self.samples[lane]
                lane_jobs.append(bcl2fastq_job)

            if final_fastq_jobs:
                lane_jobs.extend(final_fastq_jobs)

            # don't depend on notification commands
            self.add_copy_job_inputs(lane_jobs, lane)

            notification_command_start = config.param('fastq_notification_start', 'notification_command', required=False)
            if notification_command_start:
                notification_command_start = notification_command_start.format(
                    output_dir=self.output_dir,
                    number_of_mismatches=self.number_of_mismatches,
                    lane_number=self.lane_number,
                    mask=self.mask[lane],
                    technology=config.param('fastq', 'technology'),
                    run_id=self.run_id
                )
                # Use the same inputs and output of fastq job to send a notification each time the fastq job run
                job = Job(
                    [input],
                    ["notificationFastqStart." + lane + ".out"],
                    command=notification_command_start,
                    name="fastq_notification_start." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
#                lane_jobs.append(job)

            notification_command_end = config.param('fastq_notification_end', 'notification_command', required=False)
            if notification_command_end:
                notification_command_end = notification_command_end.format(
                    output_dir=self.output_dir,
                    lane_number=self.lane_number,
                    technology=config.param('fastq', 'technology'),
                    run_id=self.run_id
                )
                job = Job(
                    fastq_outputs,
                    ["notificationFastqEnd." + lane + ".out"], 
                    command=notification_command_end,
                    name="fastq_notification_end." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
#                lane_jobs.append(job)

            self.add_to_report_hash("fastq", lane_jobs)
            jobs.extend(lane_jobs)

        return jobs

    def sample_tag(self):
        """
        """
        jobs = []

        for lane in self.lanes:
            lane_jobs = []
            for readset in self.readsets[lane]:
                ouput_directory = os.path.join(self.output_dir, "Unaligned." + readset.lane, 'Project_' + readset.project_id, 'Sample_' + readset.name, "sample_tag.R1")

                lane_jobs.append(
                    concat_jobs([
                        bash.mkdir(ouput_directory),
                            tools.sh_sample_tag_summary(
                            readset.fastq1,
                            ouput_directory
                        ),
                        tools.sh_sample_tag_stats(
                            readset.fastq1,
                            ouput_directory,
                            readset.sample_tag,
                            readset.name
                    )],
                    name="sample_tag." + readset.name + "." + self.run_id + "." + str(readset.lane),
                    samples=[readset.sample]
                ))
                self.report_inputs[lane]['sample_tag'][readset.name] = os.path.join(ouput_directory, readset.name + ".sample_tag_stats.csv")

            self.add_to_report_hash("sample_tag", lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)

        return self.throttle_jobs(jobs)

    def qc_graphs(self):
        """ 
        Generate some QC Graphics and a summary XML file for each sample using 
        [BVATools](https://bitbucket.org/mugqic/bvatools/).

        Files are created in a 'qc' subfolder of the fastq directory. Examples of
        output graphic:

        - Per cycle qualities, sequence content and sequence length;
        - Known sequences (adaptors);
        - Abundant Duplicates;
        """
        jobs = []

        for lane in self.lanes:
            lane_jobs = []

            for readset in self.readsets[lane]:
                output_dir = os.path.join(os.path.dirname(readset.fastq1), "qc")
                region_name = readset.name + "_" + readset.sample_number + "_L00" + readset.lane

                file1 = readset.fastq1
                file2 = readset.fastq2
                type = "FASTQ"

                lane_jobs.append(
                    concat_jobs([
                        bash.mkdir(output_dir),
                        bvatools.readsqc(
                            file1,
                            file2,
                            type,
                            region_name,
                            output_dir
                    )],
                    name="qc." + readset.name + ".qc." + self.run_id + "." + lane,
                    samples=[readset.sample]
                ))
                self.report_inputs[lane]['qc'][readset.name] = os.path.join(output_dir, "mpsQC_" + region_name + "_stats.xml")

            self.add_to_report_hash("qc_graphs", lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)

        return self.throttle_jobs(jobs)

    def fastqc(self):
        """
        """
        jobs = []

        for lane in self.lanes:
            lane_jobs = []

            for readset in self.readsets[lane]:
                input_dict = {
                    readset.fastq1 : [
                        os.path.join(os.path.dirname(readset.fastq1), "fastqc.R1", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.fastq1))),
                        os.path.join(os.path.dirname(readset.fastq1), "fastqc.R1", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.fastq1))),
                    ],
                    readset.index_fastq1 : [
                        os.path.join(os.path.dirname(readset.index_fastq1), "fastqc.I1", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.index_fastq1))),
                        os.path.join(os.path.dirname(readset.index_fastq1), "fastqc.I1", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.index_fastq1))),
                    ]
                }
                if readset.run_type == "PAIRED_END":
                    input_dict[readset.fastq2] = [
                        os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.fastq2))),
                        os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.fastq2))),
                    ]
                    if readset.index_type == "DUALINDEX":
                        input_dict[readset.index_fastq2] = [
                            os.path.join(os.path.dirname(readset.index_fastq2), "fastqc.I2", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.index_fastq2))),
                            os.path.join(os.path.dirname(readset.index_fastq2), "fastqc.I2", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.index_fastq2))),
                        ]
                for input, outputs in input_dict.items():
                    job_suffix = re.search('.*fastqc\.([IR][12])', os.path.dirname(outputs[0])).group(1)
                    lane_jobs.append(
                        concat_jobs([
                            bash.mkdir(
                                os.path.dirname(outputs[0]),
                                remove=True
                            ),
                            fastqc.fastqc(
                                input,
                                None,
                                outputs,
                                adapter_file=None,
                                use_tmp=True
                        )],
                        name="fastqc." + readset.name + "_" + job_suffix,
                        samples=[readset.sample]
                    ))

            self.add_to_report_hash("fastqc", lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)
        return self.throttle_jobs(jobs)

    def blast(self):
        """ 
        Run blast on a subsample of the reads of each sample to find the 20 most
        frequent hits.

        The `runBlast.sh` tool from MUGQIC Tools is used. The number of reads to
        subsample can be configured by sample or for the whole lane. The output will be
        in the `Blast_sample` folder, under the Unaligned folder.
        """
        jobs = []

        nb_blast_to_do = config.param('blast', 'nb_blast_to_do', type="posint")
        is_nb_blast_per_lane = config.param('blast', 'is_nb_for_whole_lane', type="boolean")

        if is_nb_blast_per_lane:
            nb_blast_to_do = int(nb_blast_to_do) // len(self.readsets)

        nb_blast_to_do = max(1, nb_blast_to_do)

        for lane in self.lanes:
            lane_jobs = []

            self.report_inputs[lane]['blast'] = {}

            for readset in self.readsets[lane]:
                output_prefix = os.path.join(
                    self.output_dir,
                    "Unaligned." + readset.lane,
                    "Blast_sample",
                    readset.name + "_" + readset.sample_number + "_L00" + readset.lane
                )
                output = output_prefix + '.R1.RDP.blastHit_20MF_species.txt'
                self.report_inputs[lane]['blast'][readset.name] = os.path.join(output)
                current_jobs = [
                    bash.mkdir(os.path.dirname(output))
                ]

                fasta_file = output_prefix + ".R1.subSampled_{nb_blast_to_do}.fasta".format(
                    nb_blast_to_do=nb_blast_to_do
                )
                result_file = output_prefix + ".R1.subSampled_{nb_blast_to_do}.blastres".format(
                    nb_blast_to_do=nb_blast_to_do
                )

                inputs = [readset.fastq1, readset.fastq2]
                command = "runBlast.sh " + str(nb_blast_to_do) + " " + output_prefix + " " + readset.fastq1 + " "
                if readset.fastq2:
                    command += readset.fastq2
                current_jobs.append(
                    Job(
                        inputs,
                        [output],
                        [
                            ["blast", "module_mugqic_tools"],
                            ["blast", "module_blast"]
                        ],
                        command=command
                    )
                )

                # rRNA estimate using silva blast db, using the same subset of reads as the "normal" blast
                rrna_db = config.param('blast', 'rrna_db', required=False)
                if readset.is_rna and rrna_db:
                    rrna_result_file = result_file + "Rrna"
                    rrna_output = output_prefix + ".R1.subSampled_{nb_blast_to_do}.rrna".format(
                        nb_blast_to_do=nb_blast_to_do
                    )
                    command = """\
blastn \\
  -query {fasta_file} \\
  -db {db} \\
  -out {result_file} \\
  -perc_identity 80 \\
  -num_descriptions 1 \\
  -num_alignments 1""".format(
                        fasta_file=fasta_file,
                        result_file=rrna_result_file,
                        db=rrna_db
                    )
                    current_jobs.append(
                        Job(
                            [],
                            [],
                            [["blast", "module_blast"]],
                            command=command
                        )
                    )

                    command = """\
echo '{db}' \\
  > {output}""".format(
                        db=rrna_db,
                        output=rrna_output
                    )
                    current_jobs.append(
                        Job(
                            [],
                            [output],
                            [],
                            command=command
                        )                    
                    )

                    command = """\
grep ">" {result_file} | \\
wc -l >> {output}""".format(
                        result_file=rrna_result_file,
                        output=rrna_output
                    )
                    current_jobs.append(
                        Job(
                            [],
                            [output],
                            [],
                            command=command
                        )
                    )

                    command = """\
grep ">" {fasta_file} | \\
wc -l >> {output}""".format(
                        fasta_file=fasta_file,
                        output=rrna_output
                    )
                    current_jobs.append(
                        Job(
                            [],
                            [output],
                            [],
                            command=command
                        )
                    )

                # merge all blast steps of the readset into one job
                job = concat_jobs(
                    current_jobs,
                    name="blast." + readset.name + ".blast." + self.run_id + "." + lane,
                    samples = [readset.sample]
                )
                lane_jobs.append(job)

            self.add_to_report_hash("blast", lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)
        return self.throttle_jobs(jobs)

    def align(self):
        """
        Align the reads from the fastq file, sort the resulting .bam and create an index
        of that .bam.

        An basic aligment is performed on a sample when the `SampleRef` field of the
        Illumina sample sheet match one of the regexp in the configuration file and the
        corresponding genome (and indexes) are installed.

        `STAR` is used as a splice-junctions aware aligner when the sample
        `library_source` is `cDNA` or contains `RNA`; otherwise `BWA_mem` is used to
        align the reads.
        """
        jobs = []

        for lane in self.lanes:
            lane_jobs = []

            for readset in [readset for readset in self.readsets[lane] if readset.bam]:
                job = readset.aligner.get_alignment_job(readset)
                job.samples = [readset.sample]
                lane_jobs.append(job)

            self.add_to_report_hash("align", lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)
        return self.throttle_jobs(jobs)

    def picard_mark_duplicates(self):
        """
        Runs Picard mark duplicates on the sorted bam file.
        """
        jobs = []

        for lane in self.lanes:
            lane_jobs = []

            for readset in [readset for readset in self.readsets[lane] if readset.bam]:
                input_file_prefix = readset.bam + '.'
                input = input_file_prefix + "bam"
                output = input_file_prefix + "dup.bam"
                metrics_file = input_file_prefix + "dup.metrics"
                self.report_inputs[lane]['mark_dup'][readset.name] = metrics_file

                job = picard.mark_duplicates(
                    [input],
                    output,
                    metrics_file
                )
                job.name = "picard_mark_duplicates." + readset.name + ".dup." + self.run_id + "." + lane
                job.samples = [readset.sample]
                lane_jobs.append(job)

            self.add_to_report_hash("picard_mark_duplicates", lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)
        return self.throttle_jobs(jobs)

    def metrics(self):
        """
        This step runs a series of multiple metrics collection jobs and the output bam
        from mark duplicates.

        - Picard CollectMultipleMetrics: A collection of picard metrics that runs at the
        same time to save on I/O.
            - CollectAlignmentSummaryMetrics,
            - CollectInsertSizeMetrics,
            - QualityScoreDistribution,
            - MeanQualityByCycle,
            - CollectBaseDistributionByCycle
        - BVATools DepthOfCoverage: Using the specified `BED Files` in the sample sheet,
        calculate the coverage of each target region.
        - Picard CalculateHsMetrics: Calculates a set of Hybrid Selection specific
        metrics from the BAM file. The bait and interval list is automatically created
        from the specicied `BED Files`.
        """
        jobs = []

        for lane in self.lanes:
            lane_jobs = []

            for readset in [readset for readset in self.readsets[lane] if readset.bam]:
                job_list = readset.aligner.get_metrics_jobs(readset)
                for job in job_list:
                    job.samples = [readset.sample]
                lane_jobs.extend(job_list)

                if readset.is_rna:
                    self.report_inputs[lane]['align'][readset.name] = [
                        readset.bam + '.metrics.alignment_summary_metrics',
                        readset.bam + '.metrics.insert_size_metrics',
                        os.path.join(os.path.dirname(readset.bam), readset.sample.name + "." + readset.library + '.rnaseqc.sorted.dup.metrics.tsv'),
                        readset.bam + '.metrics.rRNA.tsv'
                    ]
                else:
                    self.report_inputs[lane]['align'][readset.name] = [
                        readset.bam + '.metrics.alignment_summary_metrics',
                        readset.bam + '.metrics.insert_size_metrics',
                        readset.bam + ".dup.metrics",
                        readset.bam + ".metrics.verifyBamId.tsv",
                        readset.bam + ".metrics.targetCoverage.txt"
                    ]

            self.add_to_report_hash("metrics", lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)

        return self.throttle_jobs(jobs)

    def md5(self):
        """
        Create md5 checksum files for the fastq, bam and bai using the system 'md5sum'
        util.

        One checksum file is created for each file.
        """
        jobs = []

        for lane in self.lanes:
            lane_jobs = []

            for readset in self.readsets[lane]:
                current_jobs = [
                    Job(
                        [readset.fastq1],
                        [readset.fastq1 + ".md5"],
                        command="md5sum -b " + readset.fastq1 + " > " + readset.fastq1 + ".md5"
                    )
                ]

                # Second read in paired-end run
                if readset.fastq2:
                    current_jobs.append(
                        Job(
                            [readset.fastq2],
                            [readset.fastq2 + ".md5"],
                            command="md5sum -b " + readset.fastq2 + " > " + readset.fastq2 + ".md5"
                        )
                    )

                # Alignment files
                if readset.bam:
                    current_jobs.append(
                        Job(
                            [readset.bam + ".bam"],
                            [readset.bam + ".bam.md5"],
                            command="md5sum -b " + readset.bam + ".bam" + " > " + readset.bam + ".bam.md5"
                        )
                    )
                    current_jobs.append(
                        Job(
                            [],
                            [readset.bam + ".bai.md5"],
                            command="md5sum -b " + readset.bam + ".bai" + " > " + readset.bam + ".bai.md5"
                        )
                    )

                job = concat_jobs(
                    current_jobs,
                    name="md5." + readset.name + ".md5." + self.run_id + "." + lane,
                    samples=[readset.sample]
                )

                lane_jobs.append(job)

            if config.param('md5', 'one_job', required=False, type="boolean"):
                job = concat_jobs(
                    lane_jobs,
                    name="md5." + self.run_id + "." + lane
                )
                self.add_copy_job_inputs([job], lane)
                jobs.extend([job])
            else:
                self.add_copy_job_inputs(lane_jobs, lane)
                jobs.extend(lane_jobs)

        return self.throttle_jobs(jobs)

    def copy(self):
        """
        Copy processed files to another place where they can be served or loaded into a
        LIMS.

        The destination folder and the command used can be set in the configuration
        file.

        An optional notification can be sent before the copy. The command used is in the configuration file.
        """
        inputs = self.copy_job_inputs
        jobs_to_concat = []

        full_destination_folder = os.path.join(
            config.param('copy', 'destination_folder', type="dirpath"),
            self.seq_category,
            self.year,
            os.path.basename(self.run_dir) + "-" + self.seqtype
        )

        for lane in self.lanes:
            inputs = self.copy_job_inputs[lane]

            # Notification
            output1 = os.path.join(self.output_dir, "notificationProcessingComplete." + lane + ".out")
            output2 = os.path.join(self.output_dir, "notificationCopyStart." + lane + ".out")

            notification_command = config.param('copy', 'notification_command', required=False)
            if notification_command:
                job = Job(
                    inputs,
                    [output1, output2],
                    name="start_copy_notification." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
                job.command = notification_command.format(
                    technology=config.param('copy', 'technology'),
                    output_dir=self.output_dir,
                    run_id=self.run_id,
                    output1=output1,
                    output2=output2,
                    lane_number=lane
                )
                jobs_to_concat.append(job)

            # Actual copy
            output = os.path.join(full_destination_folder, "copyCompleted." + lane + ".out")

            exclude_bam = config.param('copy', 'exclude_bam', required=False, type='boolean')
            exclude_fastq_with_bam = config.param('copy', 'exclude_fastq_with_bam', required=False, type='boolean')
            if exclude_bam and exclude_fastq_with_bam:
                log.warn("Excluding both BAM and fastq files")

            excluded_files = []

            if exclude_bam or exclude_fastq_with_bam:
                for readset in [readset for readset in self.readsets[lane] if readset.bam]:
                    if exclude_bam:
                        excluded_files.append(readset.bam + ".bam*")
                        excluded_files.append(readset.bam + ".bai*")
                    if exclude_fastq_with_bam and not exclude_bam:
                        excluded_files.append(readset.fastq1)
                        if readset.fastq2:
                            excluded_files.append(readset.fastq2)

            if self.run_dir != self.output_dir:
                copy_command_run_folder = config.param('copy', 'copy_command', required=False).format(
                    exclusion_clauses="",
                    lane_number=lane,
                    run_id=self.run_id,
                    source=self.run_dir,
                    run_name=os.path.basename(self.run_dir)
                )
                jobs_to_concat.append(
                    Job(
                        inputs,
                        [output],
                        command=copy_command_run_folder,
                        samples=self.samples[lane]
                    )
                )

            copy_command_output_folder = config.param('copy', 'copy_command', required=False).format(
                exclusion_clauses="\\\n".join(
                    [" --exclude '" + excludedfile.replace(self.output_dir + os.sep, "") + "'" for excludedfile in excluded_files]),
                lane_number=self.lane_number,
                run_id=self.run_id,
                source=self.output_dir,
                run_name=os.path.basename(self.run_dir)
            )
            jobs_to_concat.append(
                Job(
                    inputs,
                    [output],
                    command=copy_command_output_folder,
                    samples=self.samples[lane]
                )
            )
            jobs_to_concat.append(
                Job(
                    command="touch " + output,
                    samples=self.samples[lane]
                )
            )

        job = concat_jobs(
            jobs_to_concat,
            name="copy." + self.run_id,
            samples=self.samples[lane]
        )
        return [job]

    def end_copy_notification(self):
        """
        Send an optional notification to notify that the copy is finished.

        The command used is in the configuration file. This step is skipped when no
        command is provided.
        """
        jobs = []

        full_destination_folder = os.path.join(
            config.param('copy', 'destination_folder', type="dirpath"),
            self.seq_category,
            self.year,
            os.path.basename(self.run_dir) + "-" + self.seqtype
        )

        for lane in self.lanes:
            input = os.path.join(full_destination_folder, "copyCompleted." + lane + ".out")
            output = os.path.join(full_destination_folder, "notificationAssociation." + lane + ".out")

            notification_command = config.param('end_copy_notification', 'notification_command', required=False)
            if notification_command:
                job = Job(
                    [input],
                    [output],
                    name="end_copy_notification." + self.run_id + "." + lane
                )
                job.command = notification_command.format(
                    technology=config.param('end_copy_notification', 'technology'),
                    output_dir=self.output_dir,
                    run_name=os.path.basename(self.run_dir),
                    run_id=self.run_id,
                    output=output,
                    lane_number=lane
                )
                job.samples = self.samples[lane]
                jobs.append(job)

        job = concat_jobs(
            jobs,
            name="end_copy_notification." + self.run_id,
            samples=self.samples[lane]
        )
        return jobs

    def report(self):
        """
        Generate en JSON file reporting the whole pipeline
        """
        jobs = []

        report_dir = os.path.join(self.output_dir, "report")

        for lane in self.lanes:
            lane_jobs = []

            sample_report_dir  = os.path.join(report_dir, "sample_json", "L00" + lane)

            for readset in self.readsets[lane]:
                output_file = os.path.join(sample_report_dir, readset.name + ".report.json")
                jobs.append(
                    concat_jobs([
                        bash.mkdir(sample_report_dir),
                        tools.run_validation_sample_report(
                            readset,
                            self.report_inputs[lane],
                            output_file
                        )],
                        name="sample_report." + readset.name + "." + self.run_id + "." + str(readset.lane)
                    )
                )

            run_validation_report_json = os.path.join(report_dir, self.run_id + "." + lane + ".run_validation_report.json")
            general_information_file = os.path.join(self.output_dir, self.run_id + "." + lane + ".general_information.json")
            with open(general_information_file, 'w') as out_json:
                json.dump(self.report_hash, out_json, indent=4)

            run_validation_inputs = []
            for sample_job in jobs:
                run_validation_inputs.extend(sample_job.output_files)

            lane_jobs.append(
                concat_jobs([
                    bash.mkdir(report_dir),
                    tools.run_validation_aggregate_report(
                        general_information_file,
                        run_validation_inputs,
                        run_validation_report_json
                    )],
                    name="report." + self.run_id + "." + lane
                )
            )

            jobs.extend(lane_jobs)

        return jobs

    #
    # Utility methods
    #

    def add_copy_job_inputs(self, jobs, lane):
        for job in jobs:
            # we first remove dependencies of the current job, since we will have a dependency on that job
            self.copy_job_inputs[lane] = [item for item in self.copy_job_inputs[lane] if item not in job.input_files]
            self.copy_job_inputs[lane].extend(job.output_files)

    def add_to_report_hash(self, step_name, jobs=[]):
        report_hash = {
            "step_name" : step_name,
            "jobs" : [{
                "job_name" : job.name,
                "input_files" : [os.path.relpath(input_file, os.path.join(self.output_dir, "report")) for input_file in job.input_files],
                "output_files" : [os.path.relpath(output_file, os.path.join(self.output_dir, "report")) for output_file in job.output_files],
                "command" : job.command,                     
                "modules" : job.modules
            } for job in jobs]
        }
        self.report_hash["steps"].append(report_hash)

    def get_sequencer_index_length(self):
        """
        Returns the total number of index cycles of the run.
        """
        return sum(index_read.nb_cycles for index_read in [read for read in self.read_infos if read.is_index])

    def get_sequencer_minimum_read_length(self):
        """
        Returns the minimum number of cycles of a real read (not indexed).
        """
        return min(read.nb_cycles for read in [read for read in self.read_infos if (not read.is_index)])

    def get_indexcycles(self):
        """
        Returns the number of cycles for each index of the run.
        """
        for read in [read for read in self.read_infos if read.is_index]:
            if read.number in [1, 2]:
                index1cycles = read.nb_cycles
            if read.number in [3, 4]:
                index2cycles = read.nb_cycles
        return [index1cycles, index2cycles]

    def get_seqtype(self):
        """
        Determine which kind of sequencing (iseq, miseq, novaseq, hiseqx, hiseq4000 or hiseq2500) was performed,
        depending on the instrument used for the run
        """

        instrument = self.instrument
        instrument_file = config.param('DEFAULT', 'instrument_list_file', type='filepath', required=False)
        if not (instrument_file and os.path.isfile(instrument_file)):
            instrument_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'instrument_list.csv')

        return subprocess.check_output("grep -m1 '"+instrument+"' %s | awk -F',' '{print $3}'" % instrument_file, shell=True).strip()

    def get_instrument(self):
        """
        Parse the RunInfo.xml file for the instrument name the run has been running on
        """
        return Xml.parse(os.path.join(self.run_dir, "RunInfo.xml")).getroot().find('Run').find('Instrument').text

    def get_experiment_name(self):
        """
        Parse the RunParameters.xml file for the experiment name of the run
        """
        return Xml.parse(os.path.join(self.run_dir, "RunParameters.xml")).getroot().find('ExperimentName').text

    def get_run_number(self):
        """
        Parse the RunInfo.xml file for the run number
        """
        return Xml.parse(os.path.join(self.run_dir, "RunInfo.xml")).getroot().find('Run').get('Number')

    def get_seq_category(self):
        if "hiseq" in self.seqtype:
            return "hiseq"
        else:
            return self.seqtype

    def validate_barcodes(self, lane):
        """
        Validate all index sequences against each other to ensure they aren't in collision according to the chosen
        number of mismatches parameter.
        """
        min_allowed_distance = (2 * self.number_of_mismatches) + 1

        validated_indexes = []
        collisions = []

        for readset in self.readsets[lane]:
            current_index = readset.index.replace('-', '')

            for candidate_index in validated_indexes:
                if distance(current_index, candidate_index) < min_allowed_distance:
                    collisions.append("'" + current_index + "' and '" + candidate_index + "'")
            validated_indexes.append(current_index)

        if len(collisions) > 0:
            _raise(SanitycheckError("Barcode collisions: " + ";".join(collisions)))

    def get_mask(self, lane):
        """
        Returns a BCL2FASTQ friendly mask of the reads cycles.

        The mask is calculated using:
            - first base and last base of index;
            - the index length in the sample sheet;
            - the number of index cycles on the sequencer;
        """
        mask = ""
        index_lengths = self.get_smallest_index_length(lane)
        index_read_count = 0
        nb_total_index_base_used = 0

        for read_info in self.read_infos:
            if len(mask) > 0:
                mask += ','
            if read_info.is_index:
                if read_info.nb_cycles >= index_lengths[index_read_count]:
                    if index_lengths[index_read_count] == 0 or self.last_index <= nb_total_index_base_used:
                        # Don't use any index bases for this read
                        mask += 'n' + str(read_info.nb_cycles)
                    else:
                        nb_n_printed = 0

                        # Ns in the beginning of the index read
                        if self.first_index > (nb_total_index_base_used + 1):
                            nb_n_printed = min(read_info.nb_cycles, self.first_index - nb_total_index_base_used - 1)
                            if nb_n_printed >= index_lengths[index_read_count]:
                                nb_n_printed = read_info.nb_cycles
                            mask += 'n' + str(nb_n_printed)

                        # Calculate the number of index bases
                        nb_index_bases_used = max(index_lengths[index_read_count] - nb_n_printed, 0)
                        nb_index_bases_used = min(self.last_index - nb_total_index_base_used - nb_n_printed, nb_index_bases_used)
                        nb_total_index_base_used += nb_index_bases_used + min(nb_n_printed, index_lengths[index_read_count])
                        if nb_index_bases_used > 0:
                            mask += 'I' + str(nb_index_bases_used)

                        # Ns at the end of the index read
                        remaining_base_count = read_info.nb_cycles - nb_index_bases_used - nb_n_printed
                        if remaining_base_count > 0:
                            mask += 'n' + str(remaining_base_count)
                index_read_count += 1
            else:
                # Normal read
                mask += 'Y' + str(read_info.nb_cycles)
        return mask

    def generate_illumina_lane_casava_sheet(self, lane):
        """
        Create a CASAVA sheet to use with the BCL2FASTQ software
        from the Clarity data.

        Only the samples of the chosen lane will be in the file.
        The sample indexes are trimmed according to the mask used.
        """

        csv_headers = [
            "FCID",
            "Lane",
            "Sample_ID",
            "Sample_Name",
            "SampleRef",
            "Index",
            "Index2",
            "Description",
            "Control",
            "Recipe",
            "Operator",
            "Sample_Project"
        ]
        csv_file = os.path.join(
            self.output_dir,
            "casavasheet." + lane + ".indexed.csv"
        )
        writer = csv.DictWriter(
            open(csv_file, 'wb'),
            delimiter=str(','),
            fieldnames=csv_headers
        )

        # add [Data] line before the actual headers
        section_header_dict = { "FCID": "[Data]" }
        writer.writerow(section_header_dict)

        writer.writeheader()

        mask = self.mask[lane]

        overmask = ""
        overindex1 = None
        overindex2 = None

        dualindex_demultiplexing = False    # This is not the sequencing demultiplexing flag, but really the index demultiplexing flag
        self._umi = False

        # barcode validation
        if re.search("I", mask):
            self.validate_barcodes(lane)

        # IDT - UMI9 in index2
        if re.search(",I17", mask):
            self._umi = True
            overmask = re.sub(",I17,", ",I8n*,", mask)
            overindex1=8
            overindex2=8

        # HaloPlex - UMI8 in index2 alone
        if sum(1 for readset in [readset for readset in self.readsets[lane] if re.search("HaloPlex", readset.protocol)]) > 0:
            if "DUALINDEX" in set([readset.index_type for readset in self.readsets[lane]]) :
                _raise(SanityCheckError("HaloPlex libraries cannot be mixed with DUAL INDEX libraries"))
            overmask=re.sub(",I8,I10,", ",I8,Y10,", mask)
            overindex1=8
            overindex2=0
            self._bcl2fastq_extra_option="--mask-short-adapter-reads 10"

        # If SINGLEINDEX only
        if "DUALINDEX" not in set([readset.index_type for readset in self.readsets[lane]]):
            if re.search("I", mask.split(",")[2]):
                split_mask = mask.split(",") if overmask == "" else overmask.split(",")
                overmask=','.join([split_mask[0], split_mask[1], "n*", split_mask[3]])
                overindex1=8
                overindex2=0

        final_mask = mask if overmask == "" else overmask
        final_index1 = self.index1cycles if overindex1 is None else overindex1
        final_index2 = self.index2cycles if overindex2 is None else overindex2

        self._mask[lane] = config.param('fastq', 'overmask') if config.param('fastq', 'overmask', required=False, type='string') else final_mask
        self._index1cycles = config.param('fastq', 'overindex1') if config.param('fastq', 'overindex1', required=False, type='int') else final_index1
        self._index2cycles = config.param('fastq', 'overindex2') if config.param('fastq', 'overindex2', required=False, type='int') else final_index2

        # If the second index exists
        if self.index2cycles != 0:
            dualindex_demultiplexing = True

        # In case of HaloPlex-like masks, R2 is actually I2 (while R3 is R2),
        # Proceed to dualindex demultiplexing
        if ''.join(i for i in self.mask if not i.isdigit()) == "Y,I,Y,Y":
            dualindex_demultiplexing = True
        output_dir = os.path.join(self.output_dir, "Unaligned." + lane)
        casava_sample_sheet = os.path.join(self.output_dir, "casavasheet." + lane + ".indexed.csv")

        count = 0
        index_per_readset = {}
        for readset in self.readsets[lane]:
            count += 1
            index_per_readset[readset.name] = readset.indexes

            for readset_index in readset.indexes:
                csv_dict = {
                    "FCID": readset.flow_cell,
                    "Lane": lane,
                    "Sample_ID": "Sample_" + readset_index['SAMPLESHEET_NAME'],
                    "Sample_Name": readset_index['SAMPLESHEET_NAME'],
                    "SampleRef": "",
                    "Index": readset_index['INDEX1'],
                    "Index2": readset_index['INDEX2'],
                    "Description": readset.description + ' - ' + readset.protocol + ' - ' + readset.library_source,
                    "Control": readset.control,
                    "Recipe": readset.recipe,
                    "Operator": readset.operator,
                    "Sample_Project": "Project_" + readset.project_id
                }
                writer.writerow(csv_dict)
        self._index_per_readset = index_per_readset

        # Create index folder if not already done
        if not os.path.exists(os.path.join(self.output_dir, "index")):
            try:
                os.makedirs(os.path.join(self.output_dir, "index"))
            except OSError as exc: # Guard against race condition
               if exc.errno != errno.EEXIST:
                   raise
        # Python3 syntax
        #os.makedirs(os.path.dirname(index_file), exist_ok=True)

        # Print to file
        index_json = os.path.join(self.output_dir, "index", self.run_id + "_" + lane + ".index_per_readset.json")
        with open(index_json, 'w') as out_json:
            out_json.write(json.dumps(index_per_readset, indent=4))

        if self.umi:
            output_dir_noindex = output_dir + ".noindex"
            casava_sample_sheet_noindex = re.sub(".indexed.", ".noindex.", casava_sample_sheet)
            writer = csv.DictWriter(
                open(casava_sample_sheet_noindex, 'wb'),
                delimiter=str(','),
                fieldnames=csv_headers
            )
            writer.writerow(
                {
                    "FCID": readset.flow_cell,
                    "Lane": self.lane_number,
                    "Sample_ID": "Sample_ALL",
                    "Sample_Name": "ALL",
                    "SampleRef": "",
                    "Index": "",
                    "Index2": "",
                    "Description": "",
                    "Control": "N",
                    "Recipe": "",
                    "Operator": "",
                    "Sample_Project": "Project_ALL"
                }
            )

    def generate_fastq_outputs(self, lane):
        bcl2fastq_outputs = []
        final_fastq_jobs = []
        count = 0
        merge = self.merge_undetermined

        output_dir = os.path.join(self.output_dir, "Unaligned." + lane)

        for readset in self.readsets[lane]:
            final_fastq_jobs_per_readset = []
            count += 1

            fastq_output_dir = os.path.join(output_dir, 'Project_' + readset.project_id, 'Sample_' + readset.name)
            # If 10X libraries : 4 indexes per sample
            if re.search("tenX", readset.library_type):
                fastq1 = os.path.basename(readset.fastq1)

                counta = count
                fastq1a=os.path.join(fastq_output_dir + "_A", re.sub(readset.name + "_S" + readset.sample_number, readset.name + "_A_S" + str(counta), fastq1))
                bcl2fastq_outputs.append(fastq1a)

                countb = count + 1
                fastq1b=os.path.join(fastq_output_dir + "_B", re.sub(readset.name + "_S" + readset.sample_number, readset.name + "_B_S" + str(countb), fastq1))
                bcl2fastq_outputs.append(fastq1b)

                countc = count + 2
                fastq1c=os.path.join(fastq_output_dir + "_C", re.sub(readset.name + "_S" + readset.sample_number, readset.name + "_C_S" + str(countc), fastq1))
                bcl2fastq_outputs.append(fastq1c)

                countd = count + 3
                fastq1d=os.path.join(fastq_output_dir + "_D", re.sub(readset.name + "_S" + readset.sample_number, readset.name + "_D_S" + str(countd), fastq1))
                bcl2fastq_outputs.append(fastq1d)

                final_fastq_jobs_per_readset.append(
                    bash.cat(
                        [
                            fastq1a,
                            fastq1b,
                            fastq1c,
                            fastq1d
                        ],
                        readset.fastq1
                ))
                # If True, then merge the 'Undetermined' reads
                if merge:
                    undet_fastq1 = os.path.join(output_dir, re.sub(readset.name, "Undetermined_S0", fastq1))
                    bcl2fastq_outputs.append(undet_fastq1)
                    final_fastq_jobs_per_readset.append(
                        bash.cat(
                            [undet_fastq1],
                            readset.fastq1,
                            append=True
                    ))

                # Add the fastq of first index
                idx_fastq1 = re.sub("_R1_", "_I1_", readset.fastq1)
                idx_fastq1a = re.sub("_R1_", "_I1_", fastq1a)
                idx_fastq1b = re.sub("_R1_", "_I1_", fastq1b)
                idx_fastq1c = re.sub("_R1_", "_I1_", fastq1c)
                idx_fastq1d = re.sub("_R1_", "_I1_", fastq1d)
                bcl2fastq_outputs.append(idx_fastq1a)
                bcl2fastq_outputs.append(idx_fastq1b)
                bcl2fastq_outputs.append(idx_fastq1c)
                bcl2fastq_outputs.append(idx_fastq1d)
                final_fastq_jobs_per_readset.append(
                    bash.cat(
                        [
                            idx_fastq1a,
                            idx_fastq1b,
                            idx_fastq1c,
                            idx_fastq1d
                        ],
                        idx_fastq1
                ))
                # If True, then merge the 'Undetermined' reads
                if merge:
                    undet_idx_fastq1 = os.path.join(output_dir, re.sub(readset.name, "Undetermined_S0", os.path.basename(idx_fastq1)))
                    bcl2fastq_outputs.append(undet_idx_fastq1)
                    final_fastq_jobs_per_readset.append(
                        bash.cat(
                            [undet_idx_fastq1],
                            idx_fastq1,
                            append=True
                    ))

                # For paired-end sequencing, do not forget the fastq of the reverse reads
                if readset.run_type == "PAIRED_END" :
                    fastq2 = os.path.basename(readset.fastq2)
                    fastq2a = os.path.join(fastq_output_dir + "_A", re.sub(readset.name + "_S" + readset.sample_number, readset.name + "_A_S" + str(counta), fastq2))
                    fastq2b = os.path.join(fastq_output_dir + "_B", re.sub(readset.name + "_S" + readset.sample_number, readset.name + "_B_S" + str(counta), fastq2))
                    fastq2c = os.path.join(fastq_output_dir + "_C", re.sub(readset.name + "_S" + readset.sample_number, readset.name + "_C_S" + str(counta), fastq2))
                    fastq2d = os.path.join(fastq_output_dir + "_D", re.sub(readset.name + "_S" + readset.sample_number, readset.name + "_D_S" + str(counta), fastq2))
                    bcl2fastq_outputs.append(fastq2a)
                    bcl2fastq_outputs.append(fastq2b)
                    bcl2fastq_outputs.append(fastq2c)
                    bcl2fastq_outputs.append(fastq2d)
                    final_fastq_jobs_per_readset.append(
                        bash.cat(
                            [
                                fastq2a,
                                fastq2b,
                                fastq2c,
                                fastq2d
                            ],
                            readset.fastq2
                    ))
                    # If True, then merge the 'Undetermined' reads
                    if merge:
                        undet_fastq2 = os.path.join(output_dir, re.sub(readset.name, "Undetermined_S0", fastq2))

                        bcl2fastq_outputs.append(undet_fastq2)
                        final_fastq_jobs_per_readset.append(
                            bash.cat(
                                [undet_fastq2],
                                readset.fastq2,
                                append=True
                        ))

                    if readset.index_type == "DUALINDEX" :
                        # For dual index demultiplexing, do not forget the fastq of the second index
                        idx_fastq2 = re.sub("_R2_", "_I2_", readset.fastq2)
                        idx_fastq2a = re.sub("_R1_", "_I1_", fastq2a)
                        idx_fastq2b = re.sub("_R1_", "_I1_", fastq2b)
                        idx_fastq2c = re.sub("_R1_", "_I1_", fastq2c)
                        idx_fastq2d = re.sub("_R1_", "_I1_", fastq2d)
                        bcl2fastq_outputs.append(idx_fastq2a)
                        bcl2fastq_outputs.append(idx_fastq2b)
                        bcl2fastq_outputs.append(idx_fastq2c)
                        bcl2fastq_outputs.append(idx_fastq2d)
                        final_fastq_jobs_per_readset.append(
                            bash.cat(
                                [
                                    idx_fastq2a,
                                    idx_fastq2b,
                                    idx_fastq2c,
                                    idx_fastq2d
                                ],
                                idx_fastq2
                        ))
                        # If True, then merge the 'Undetermined' reads
                        if merge:
                            undet_idx_fastq2 = os.path.join(output_dir, re.sub(readset.name, "Undetermined_S0", os.path.basename(idx_fastq2)))
                            bcl2fastq_outputs.append(re.sub(undet_idx_fastq2))
                            final_fastq_jobs_per_readset.append(
                                bash.cat(
                                    [undet_idx_fastq2],
                                    idx_fastq2,
                                    append=True
                            ))

                # For HaloPlex-like masks, R3 fastq is created
                if ''.join(i for i in self.mask if not i.isdigit()) == "Y,I,Y,Y":
                    fastq3 = re.sub("_R1_", "_R3_", os.path.basename(readset.fastq1))
                    fastq3a = os.path.join(fastq_output_dir + "_A", re.sub(readset.name + "_S" + readset.sample_number, readset.name + "_A_S" + str(counta), fastq3))
                    fastq3b = os.path.join(fastq_output_dir + "_B", re.sub(readset.name + "_S" + readset.sample_number, readset.name + "_B_S" + str(counta), fastq3))
                    fastq3c = os.path.join(fastq_output_dir + "_C", re.sub(readset.name + "_S" + readset.sample_number, readset.name + "_C_S" + str(counta), fastq3))
                    fastq3d = os.path.join(fastq_output_dir + "_D", re.sub(readset.name + "_S" + readset.sample_number, readset.name + "_D_S" + str(counta), fastq3))
                    bcl2fastq_outputs.append(fastq3a)
                    bcl2fastq_outputs.append(fastq3b)
                    bcl2fastq_outputs.append(fastq3c)
                    bcl2fastq_outputs.append(fastq3d)
                    final_fastq_jobs_per_readset.append(
                        bash.cat(
                            [
                                fastq3a,
                                fastq3b,
                                fastq3c,
                                fastq3d
                            ],
                            re.sub("_R1_", "_R3_", readset.fastq1)
                    ))
                    # If True, then merge the 'Undetermined' reads
                    if merge:
                        undet_fastq3 = os.path.join(output_dir, re.sub(readset.name, "Undetermined_S0", fastq3))
                        bcl2fastq_outputs.append(undet_fastq3)
                        final_fastq_jobs_per_readset.append(
                            bash.cat(
                                [undet_fastq3],
                                re.sub("_R1_", "_R3_", readset.fastq1),
                                append=True
                        ))

                count = countd

            # not a 10X library : 1 index per sample
            else:
                bcl2fastq_outputs.append(readset.fastq1)
                # If ask to merge the Undeternined reads
                if merge:
                    undet_fastq1 = os.path.join(output_dir, re.sub(readset.name + "_S" + readset.sample_number, "Undetermined_S0", os.path.basename(readset.fastq1)))
                    bcl2fastq_outputs.append(undet_fastq1)
                    final_fastq_jobs_per_readset.append(
                        bash.cat(
                            [undet_fastq1],
                            readset.fastq1,
                            append=True
                    ))

                idx_fastq1 = re.sub("_R1_", "_I1_", readset.fastq1)
                bcl2fastq_outputs.append(idx_fastq1)
                # If ask to merge the Undeternined reads
                if merge:
                    undet_idx_fastq1 = os.path.join(output_dir, re.sub(readset.name + "_S" + readset.sample_number, "Undetermined_S0", os.path.basename(idx_fastq1)))
                    bcl2fastq_outputs.append(undet_idx_fastq1)
                    final_fastq_jobs_per_readset.append(
                        bash.cat(
                            [undet_idx_fastq1],
                            idx_fastq1,
                            append=True
                    ))

                # For paired-end sequencing, do not forget the fastq of the reverse reads
                if readset.run_type == "PAIRED_END" :
                    bcl2fastq_outputs.append(readset.fastq2)
                    # If True, then merge the 'Undetermined' reads
                    if merge:
                        undet_fastq2 = os.path.join(output_dir, re.sub(readset.name + "_S" + readset.sample_number, "Undetermined_S0", os.path.basename(readset.fastq2)))
                        bcl2fastq_outputs.append(undet_fastq2)
                        final_fastq_jobs_per_readset.append(
                            bash.cat(
                                [undet_fastq2],
                                readset.fastq2,
                                append=True
                        ))

                    if readset.index_type == "DUALINDEX" :
                        # For dual index multiplexing, do not forget the fastq of the second index
                        idx_fastq2 = re.sub("_R2_", "_I2_", readset.fastq2)
                        bcl2fastq_outputs.append(idx_fastq2)
                        # If True, then merge the 'Undetermined' reads
                        if merge:
                            undet_idx_fastq2 = os.path.join(output_dir, re.sub(readset.name + "_S" + readset.sample_number, "Undetermined_S0", os.path.basename(idx_fastq2)))
                            bcl2fastq_outputs.append(undet_idx_fastq2)
                            final_fastq_jobs_per_readset.append(
                                bash.cat(
                                    [undet_idx_fastq2],
                                    idx_fastq2,
                                    append=True
                            ))

                # For HaloPlex-like masks, do the necessary name swapping
                if ''.join(i for i in self.mask if not i.isdigit()) == "Y,I,Y,Y":
                    fastq3 = re.sub("_R1_", "_R3_", readset.fastq1)
                    bcl2fastq_outputs.append(fastq3)
                    # If True, then merge the 'Undetermined' reads
                    if merge:
                        undet_fastq3 = os.path.join(output_dir, re.sub(readset.name+ "_S" + readset.sample_number, "Undetermined_S0", os.path.basename(fastq3)))
                        bcl2fastq_outputs.append(undet_fastq3)
                        final_fastq_jobs_per_readset.append(
                            bash.cat(
                                [undet_fastq3],
                                fastq3,
                                append=True
                        ))

            # For HaloPlex-like masks, do the necessary name swapping
            if ''.join(i for i in self.mask if not i.isdigit()) == "Y,I,Y,Y":
                final_fastq_jobs_per_readset.append(
                    concat_jobs([
                        bash.mv(
                            fastq2,
                            idx_fastq2
                        ),
                        bash.mv(
                            fastq3,
                            fastq2
                    )]
                ))

            if len(final_fastq_jobs_per_readset) != 0:
                final_fastq_jobs.append(
                    concat_jobs(
                        final_fastq_jobs_per_readset,
                        name="build_final_fastq."+readset.name,
                        samples=[readset.sample]
                    )
                )

        return bcl2fastq_outputs, final_fastq_jobs

    def has_single_index(self, lane):
        """ 
        Returns True when there is at least one sample on the lane that doesn't use double-indexing or we only have
        one read of indexes.
        """
        return len([readset for readset in self.readsets[lane] if ("-" not in readset.index)]) > 0 or\
               len([read for read in self.read_infos[lane] if read.is_index]) < 2

    def get_smallest_index_length(self, lane):
        """
        Returns a list (for each index read of the lane) of the minimum between the number of index cycle on the
        sequencer and all the index lengths.
        """
        run_index_lengths = [r.nb_cycles for r in self.read_infos if r.is_index] # from RunInfo

        if len(run_index_lengths) == 0 and len(self.readsets[lane]) > 1:
            _raise(SanitycheckError("Multiple samples on lane '" + lane + "', but no indexes were read from the sequencer."))

        # loop on all index reads, to compare with samples index length
        for i in range(len(run_index_lengths)):
            min_sample_index_length = 0
            try:
                min_sample_index_length = min(len(readset.index.split("-")[i])
                                              for readset in self.readsets[lane]
                                              if (len(readset.index.split("-")) > i and len(readset.index.split("-")[i]) > 0))
            except ValueError:
                pass  # we don't have a sample with this Ith index read, use the 0 already set

            empty_index_list = [ readset for readset in self.readsets[lane] if (len(readset.index.split("-")) <= i or len(readset.index.split("-")[i]) == 0) ]
        
            if len(empty_index_list):
                # we have samples without this Ith index read, so we skip it
                min_sample_index_length = 0

            run_index_lengths[i] = min(min_sample_index_length, run_index_lengths[i])

        return run_index_lengths

    def parse_run_info_file(self):
        """
        Parse the RunInfo.xml file of the run and returns the list of RunInfoRead objects
        """
        reads = Xml.parse(os.path.join(self.run_dir, "RunInfo.xml")).getroot().find('Run').find('Reads')
        return [ RunInfoRead(int(r.get("Number")), int(r.get("NumCycles")), r.get("IsIndexedRead") == "Y") for r in reads.iter('Read') ]

    def load_readsets(self, lane):
        """
        Parse the sample sheet and return a list of readsets.
        """

        return parse_illumina_raw_readset_files(
            self.output_dir,
            self.run_dir,
            "PAIRED_END" if self.is_paired_end else "SINGLE_END",
            self.readset_file,
            lane,
            config.param('DEFAULT', 'genomes_home', type="dirpath"),
            self.get_sequencer_minimum_read_length(),
            self.index1cycles,
            self.index2cycles,
            self.seqtype
        )

    def submit_jobs(self):
        super(IlluminaRunProcessing, self).submit_jobs()

    def throttle_jobs(self, jobs):
        """
        Group jobs of the same task (same name prefix) if they exceed the configured threshold number.
        """

        max_jobs_per_step = config.param('default', 'max_jobs_per_step', required=False, type="int")
        jobs_by_name = OrderedDict()
        reply = []

        # group jobs by task (name)
        for job in jobs:
            jobs_by_name.setdefault(job.name.split(".", 1)[0], []).append(job)

        # loop on all task
        for job_name in jobs_by_name:
            current_jobs = jobs_by_name[job_name]
            if max_jobs_per_step and 0 < max_jobs_per_step < len(current_jobs):
                # we exceed the threshold, we group using 'number_task_by_job' jobs per group
                number_task_by_job = int(math.ceil(len(current_jobs) / max_jobs_per_step))
                merged_jobs = []
                for x in range(max_jobs_per_step):
                    if x * number_task_by_job < len(current_jobs):
                        merged_jobs.append(
                            concat_jobs(
                                current_jobs[x * number_task_by_job:min((x + 1) * number_task_by_job, len(current_jobs))],
                                name=job_name + "." + str(x + 1) + "." + self.run_id
                            )
                        )
                reply.extend(merged_jobs)
            else:
                reply.extend(current_jobs)
        return reply

    @property
    def steps(self):
        return [
            self.index,
            self.fastq,
            self.sample_tag,
            self.qc_graphs,
            self.fastqc,
            self.blast,
            self.align,
            self.picard_mark_duplicates,
            self.metrics,
            self.md5,
            self.report,
            self.copy,
            self.end_copy_notification
        ]

def distance(
    str1,
    str2
    ):
    """
    Returns the hamming distance. http://code.activestate.com/recipes/499304-hamming-distance/#c2
    """
    return sum(itertools.imap(unicode.__ne__, str1, str2))

if __name__ == '__main__':

    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        IlluminaRunProcessing()

