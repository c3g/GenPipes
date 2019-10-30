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
import os
import re
import sys
import itertools
import subprocess
import xml.etree.ElementTree as Xml
import math
import subprocess
import csv

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from bfx.readset import parse_illumina_raw_readset_files

from bfx import bvatools
from bfx import picard
from bfx import tools
from bfx import run_processing_tools
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
        self.copy_job_inputs = []
        self.argparser.add_argument("-d", "--run", help="Run directory", required=False, dest="run_dir")
        self.argparser.add_argument("--lane", help="Lane number", type=int, required=False, dest="lane_number")
        self.argparser.add_argument("-r", "--readsets", help="If type='clarity', clarity event file (mandatory).\nIf type='nanuq', nanuq readset file. The default file is 'run.nanuq.csv' in the output folder. Will be automatically downloaded if not present.", type=file, required=False)
        self.argparser.add_argument("-i", help="Illumina casava sheet.\nIf type='nanuq', the default file is 'SampleSheet.nanuq.csv' in the output folder. Will be automatically downloaded if not present.\nIgnored if type='clarity'", type=file, required=False, dest="casava_sheet_file")
        self.argparser.add_argument("-x", help="First index base to use for demultiplexing (inclusive). The index from the sample sheet will be adjusted according to that value.", type=int, required=False, dest="first_index")
        self.argparser.add_argument("-y", help="Last index base to use for demultiplexing (inclusive)", type=int, required=False, dest="last_index")
        self.argparser.add_argument("-m", help="Number of index mistmaches allowed for demultiplexing (default 1). Barcode collisions are always checked.", type=int, required=False, dest="number_of_mismatches")
        self.argparser.add_argument("-w", "--force-download", help="If type='nanuq', force the download of the samples sheets (default: false)\nIgnored if type='clarity'", action="store_true", dest="force_download")
        self.argparser.add_argument("-t", "--type", help="Used LIMS to produce the run (default: clarity)", choices = ["clarity", "nanuq"], default="clarity")

        super(IlluminaRunProcessing, self).__init__(protocol)

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = self.load_readsets()
            self.generate_illumina_lane_sample_sheet
        return self._readsets

    @property
    def is_paired_end(self):
        if not hasattr(self, "_is_paired_end"):
            self._is_paired_end = len([read_info for read_info in self.read_infos if (not read_info.is_index)]) > 1
        return self._is_paired_end

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
            return self.args.run_dir
        else:
            _raise(SanitycheckError("Error: missing '-d/--run' option!"))

    @property
    def lane_number(self):
        if self.args.lane_number:
            return self.args.lane_number
        else:
            _raise(SanitycheckError("Error: missing '--lane' option!"))

    @property
    def casava_sheet_file(self):
        return self.args.casava_sheet_file.name if self.args.casava_sheet_file else os.path.join(self.output_dir, "SampleSheet.nanuq.csv")

    @property
    def nanuq_readset_file(self):
        return self.args.readsets.name if self.args.readsets else os.path.join(self.output_dir, "run.nanuq.csv")

    @property
    def clarity_event_file(self):
        if self.args.readsets:
            return self.args.readsets.name
        elif self.protocol == "clarity":
            _raise(SanitycheckError("Error: missing '-r/--readsets' argument !"))

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
    def mask(self):
        if not hasattr(self, "_mask"):
            self._mask = self.get_mask()
        return self._mask

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
    def indexes_from_lims(self):
        # Define in generate_clarity_sample_sheet() 
        return self._indexes_from_lims

    @property
    def seqtype(self):
        if not hasattr(self, "_seqtype"):
            self._seqtype = self.get_seqtype()
        return self._seqtype

    @property
    def instrument(self):
        if not hasattr(self, "_instrument"):
            self._instrument = self.get_instrument()
        return self._instrument

    @property
    def read_infos(self):
        if not hasattr(self, "_read_infos"):
            self._read_infos = self.parse_run_info_file()
        return self._read_infos

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
            input = os.path.join(self.run_dir, "RunInfo.xml")
            output = os.path.join(self.output_dir, "index",  self.run_id + "_" + str(self.lane_number) + '.metrics')
            basecalls_dir = os.path.join(self.output_dir, "index", "BaseCalls")

            jobs.append(
                concat_jobs([
                    bash.mkdir(basecalls_dir),
                    bash.ln(
                        os.path.join(self.run_dir, "Data", "Intensities", "BaseCalls", "L00" + self.lane_number),
                        os.path.join(basecalls_dir, "L00" + self.lane_number)
                    ) if not os.path.islink(os.path.join(basecalls_dir, "L00" + lane)) else None,
                    bash.ln(
                        os.path.join(self.run_dir, "Data", "Intensities", "s.locs"),
                        os.path.join(basecalls_dir, "..", "index", "s.locs")
                    ) if not os.path.islink(os.path.join(basecalls_dir, "..", "index", "s.locs")) else  None
                    run_processing_tools.index(
                        input,
                        output,
                        basecalls_dir,
                        self.number_of_mismatches,
                        self.lane_number,
                        mask
                )],
                name="index." + self.run_id + "." + str(self.lane_number)
                samples=self.samples
            ))

        # Prepare Index Validation jobs
        indexes_from_lims = self.lims_index

        self.add_copy_job_inputs(jobs)
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

        if self.protocol == "nanuq":
            jobs = self.nanuq_fastq()
        else:
            jobs = self.clarity_fastq()
        
        return jobs

    def nanuq_fastq(self):
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

        input = self.bcl2fastq_job_input

        [
            fastq_outputs,
            output_dir,
            bcl2fastq_extra_option
            
        ] = self.prepare_bcl2fastq_inputs()

        bcl2fastq_job = run_processing_tools.bcl2fastq(
            input,
            fastq_outputs,
            output_dir,
            casava_sample_sheet,
            self.run_dir,
            self.lane_number,
            bcl2fastq_extra_option,
            demultiplex=True,
            mismatches=self.number_of_mismatches,
            mask=mask,
        )

        

        aggregate_fastq_jobs = None
        for readset in self.readsets:
            if re.search("tenX", readset.library_type) or merge_undetermined:
                aggregate_fastq_jobs = concat_jobs([
                    aggregate_fastq_jobs,
                    run_processing_tools.aggregate_fastqs(
                        readset,
                        merge_undetermined,
                        mask
                ])

        if generate_umi:

            

            output_dir_noindex = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number) + ".noindex")
            casava_sample_sheet_noindex = os.path.join(self.output_dir, "casavasheet." + str(self.lane_number) + ".noindex.csv")

            

            jobs.append(
                concat_jobs([
                    bcl2fastq_job,
                    run_processing_tools.aggregate_fastqs(
                        merge_undetermined
                )]
            ))
            
                    run_processing_tools.bcl2fastq(
                        input,
                        fastq_outputs,
                        output_dir_noindex,
                        casava_sample_sheet_noindex,
                        self.run_dir,
                        self.lane_number,
                        bcl2fastq_extra_option
                    ),
                name="fastq." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            ))
        
        else:
            output_dir = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number))
            casava_sheet_prefix = config.param('fastq', 'casava_sample_sheet_prefix')
            casava_sample_sheet = os.path.join(self.output_dir, "casavasheet." + str(self.lane_number) + ".indexed.csv")

            jobs.append(
                run_processing_tools.bcl2fastq(
                    input,
                    fastq_outputs,
                    output_dir,
                    sample_sheet,
                    self.run_dir,
                    self.lane_number,
                    bcl2fastq_extra_option,
                    demultiplex=True,
                    mismatches=self.number_of_mismatches,
                    mask=mask,
            ))
                name="fastq." + self.run_id + "." + str(self.lane_number),
        samples=self.samples
            


        # don't depend on notification commands
        self.add_copy_job_inputs(jobs)

        notification_command_start = config.param('fastq_notification_start', 'notification_command', required=False)
        if notification_command_start:
            notification_command_start = notification_command_start.format(
                output_dir=self.output_dir,
                number_of_mismatches=self.number_of_mismatches if demultiplexing else "-",
                lane_number=self.lane_number,
                mask=mask if demultiplexing else "-",
                technology=config.param('fastq', 'technology'),
                run_id=self.run_id
            )
            # Use the same inputs and output of fastq job to send a notification each time the fastq job run
            job = Job(
                [input],
                ["notificationFastqStart." + str(self.lane_number) + ".out"],
                command=notification_command_start,
                name="fastq_notification_start." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            )
            jobs.append(job)

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
                ["notificationFastqEnd." + str(self.lane_number) + ".out"], 
                command=notification_command_end,
                name="fastq_notification_end." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            )
            jobs.append(job)

        return jobs

    def clarity_fastq(self):
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

        input = self.clarity_event_file

        fastq_outputs = self.generate_fastq_outputs(merge)

        bcl2fastq_job = run_processing_tools.bcl2fastq(
            input,
            fastq_outputs,
            output_dir,
            casava_sample_sheet,
            self.run_dir,
            self.lane_number,
            bcl2fastq_extra_option,
            demultiplex=True,
            mismatches=self.number_of_mismatches,
            mask=mask,
        )

        aggregate_fastq_jobs = None
        for readset in self.readsets:
            if re.search("tenX", readset.library_type) or merge_undetermined:
                aggregate_fastq_jobs = concat_jobs([
                    aggregate_fastq_jobs,
                    run_processing_tools.aggregate_fastqs(
                        readset,
                        merge_undetermined,
                        mask
                ])

        if generate_umi:

            output_dir_noindex = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number) + ".noindex")
            casava_sample_sheet_noindex = os.path.join(self.output_dir, "casavasheet." + str(self.lane_number) + ".noindex.csv")

            jobs.append(
                concat_jobs([
                    bcl2fastq_job,
                    run_processing_tools.aggregate_fastqs(
                        merge_undetermined
                )]
            ))

            run_processing_tools.bcl2fastq(
                input,
                fastq_outputs,
                output_dir_noindex,
                casava_sample_sheet_noindex,
                self.run_dir,
                self.lane_number,
                bcl2fastq_extra_option
            ),
            name="fastq." + self.run_id + "." + str(self.lane_number),
            samples=self.samples

        else:
            output_dir = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number))
            casava_sheet_prefix = config.param('fastq', 'casava_sample_sheet_prefix')
            casava_sample_sheet = os.path.join(self.output_dir, "casavasheet." + str(self.lane_number) + ".indexed.csv")

            jobs.append(
                run_processing_tools.bcl2fastq(
                    input,
                    fastq_outputs,
                    output_dir,
                    sample_sheet,
                    self.run_dir,
                    self.lane_number,
                    bcl2fastq_extra_option,
                    demultiplex=True,
                    mismatches=self.number_of_mismatches,
                    mask=mask,
            ))
            name="fastq." + self.run_id + "." + str(self.lane_number),
            samples=self.samples

        # don't depend on notification commands
        self.add_copy_job_inputs(jobs)

        notification_command_start = config.param('fastq_notification_start', 'notification_command', required=False)
        if notification_command_start:
            notification_command_start = notification_command_start.format(
                output_dir=self.output_dir,
                number_of_mismatches=self.number_of_mismatches if demultiplexing else "-",
                lane_number=self.lane_number,
                mask=mask if demultiplexing else "-",
                technology=config.param('fastq', 'technology'),
                run_id=self.run_id
            )
            # Use the same inputs and output of fastq job to send a notification each time the fastq job run
            job = Job(
                [input],
                ["notificationFastqStart." + str(self.lane_number) + ".out"],
                command=notification_command_start,
                name="fastq_notification_start." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            )
            jobs.append(job)

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
                ["notificationFastqEnd." + str(self.lane_number) + ".out"], 
                command=notification_command_end,
                name="fastq_notification_end." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            )
            jobs.append(job)

        return jobs

    def clarity_fastq(self):
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
        

        input = self.bcl2fastq_job_input

        [
            fastq_outputs,
            output_dir,
            bcl2fastq_extra_option
            
        ] = self.prepare_bcl2fastq_inputs()

        bcl2fastq_job = run_processing_tools.bcl2fastq(
            input,
            fastq_outputs,
            output_dir,
            casava_sample_sheet,
            self.run_dir,
            self.lane_number,
            bcl2fastq_extra_option,
            demultiplex=True,
            mismatches=self.number_of_mismatches,
            mask=mask,
        )

        

        aggregate_fastq_jobs = None
        for readset in self.readsets:
            if re.search("tenX", readset.library_type) or merge_undetermined:
                aggregate_fastq_jobs = concat_jobs([
                    aggregate_fastq_jobs,
                    run_processing_tools.aggregate_fastqs(
                        readset,
                        merge_undetermined,
                        mask
                ])

        if generate_umi:

            

            output_dir_noindex = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number) + ".noindex")
            casava_sample_sheet_noindex = os.path.join(self.output_dir, "casavasheet." + str(self.lane_number) + ".noindex.csv")

            

            jobs.append(
                concat_jobs([
                    bcl2fastq_job,
                    run_processing_tools.aggregate_fastqs(
                        merge_undetermined
                )]
            ))
            
                    run_processing_tools.bcl2fastq(
                        input,
                        fastq_outputs,
                        output_dir_noindex,
                        casava_sample_sheet_noindex,
                        self.run_dir,
                        self.lane_number,
                        bcl2fastq_extra_option
                    ),
                name="fastq." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            ))
        
        else:
            output_dir = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number))
            casava_sheet_prefix = config.param('fastq', 'casava_sample_sheet_prefix')
            casava_sample_sheet = os.path.join(self.output_dir, "casavasheet." + str(self.lane_number) + ".indexed.csv")

            jobs.append(
                run_processing_tools.bcl2fastq(
                    input,
                    fastq_outputs,
                    output_dir,
                    sample_sheet,
                    self.run_dir,
                    self.lane_number,
                    bcl2fastq_extra_option,
                    demultiplex=True,
                    mismatches=self.number_of_mismatches,
                    mask=mask,
            ))
                name="fastq." + self.run_id + "." + str(self.lane_number),
        samples=self.samples
            


        # don't depend on notification commands
        self.add_copy_job_inputs(jobs)

        notification_command_start = config.param('fastq_notification_start', 'notification_command', required=False)
        if notification_command_start:
            notification_command_start = notification_command_start.format(
                output_dir=self.output_dir,
                number_of_mismatches=self.number_of_mismatches if demultiplexing else "-",
                lane_number=self.lane_number,
                mask=mask if demultiplexing else "-",
                technology=config.param('fastq', 'technology'),
                run_id=self.run_id
            )
            # Use the same inputs and output of fastq job to send a notification each time the fastq job run
            job = Job(
                [input],
                ["notificationFastqStart." + str(self.lane_number) + ".out"],
                command=notification_command_start,
                name="fastq_notification_start." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            )
            jobs.append(job)

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
                ["notificationFastqEnd." + str(self.lane_number) + ".out"], 
                command=notification_command_end,
                name="fastq_notification_end." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            )
            jobs.append(job)

        return jobs

    def sample_tag(self):
        """
        """
        job = []

        for readset in self.readsets:
            ouput_directory = os.path.join(self.output_dir, "Unaligned." + readset.lane, 'Project_' + readset.project, 'Sample_' + readset.name, "sample_tag.R1")                                           
            jobs.append(
                concat_jobs([
                    bash.mkdir(ouput_directory),
                    tools.sh_sample_tag_summary(
                        readset.fastq1,
                        ouput_directory
                    ),
                    tools.sh_sample_tag_stats(
                        os.path.join(ouput_directory, fastq1)
                        ouput_directory,
                        expected_sample_tag,
                        readset.name
                )],
                name="sample_tag." + readset.name + "." + self.run_id + "." + str(readset.lane)
                samples=[readset.sample]
            ))

        self.add_copy_job_inputs(jobs)
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
        for readset in [readset for readset in self.readsets if readset.bam]:
            job = readset.aligner.get_alignment_job(readset)
            job.samples = [readset.sample]
            jobs.append(job)

        self.add_copy_job_inputs(jobs)
        return self.throttle_jobs(jobs)

    def picard_mark_duplicates(self):
        """
        Runs Picard mark duplicates on the sorted bam file.
        """
        jobs = []
        for readset in [readset for readset in self.readsets if readset.bam]:
            input_file_prefix = readset.bam + '.'
            input = input_file_prefix + "bam"
            output = input_file_prefix + "dup.bam"
            metrics_file = readset.bam + ".dup.metrics"

            job = picard.mark_duplicates(
                [input],
                output,
                metrics_file
            )
            job.name = "picard_mark_duplicates." + readset.name + ".dup." + self.run_id + "." + str(self.lane_number)
            job.samples = [readset.sample]
            jobs.append(job)

        self.add_copy_job_inputs(jobs)
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
        for readset in [readset for readset in self.readsets if readset.bam]:
            job_list = readset.aligner.get_metrics_jobs(readset)
            for job in job_list:
                job.samples = [readset.sample]
            jobs.extend(job_list)
        self.add_copy_job_inputs(jobs)
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

        for readset in self.readsets:
            output_prefix = os.path.join(
                self.output_dir,
                "Unaligned." + readset.lane,
                "Blast_sample",
                readset.name + "_" + readset.sample_number + "_L00" + readset.lane
            )
            output = output_prefix + '.R1.RDP.blastHit_20MF_species.txt'
            current_jobs = [
                bash.mkdir(os.path.dirname(output))
            ]

            fasta_file = output_prefix + ".R1.subSampled_{nb_blast_to_do}.fasta".format(
                nb_blast_to_do=nb_blast_to_do
            )
            result_file = output_prefix + ".R1.subSampled_{nb_blast_to_do}.blastres".format(
                nb_blast_to_do=nb_blast_to_do
            )

            if readset.bam:
                input = readset.bam + ".bam"

                # count the read that aren't marked as secondary alignment and calculate the ratio of reads to subsample
                command = """\
subsampling=$(samtools view -F 0x0180 {input} | \\
wc -l | \\
awk -v nbReads={nb_blast_to_do} \\
  '{{x=sprintf("%.4f", nbReads/$1); if (x == "0.0000") print "0.0001"; else print x}}')""".format(
                    input=input,
                    nb_blast_to_do=nb_blast_to_do
                )
                current_jobs.append(
                    Job(
                        [input],
                        [],
                        [["blast", "module_samtools"]],
                        command=command
                    )
                )

                # subsample the reads and output to a temp fasta
                command = """samtools view \\
  -s $subsampling \\
  -F 0x0180 \\
  {input} | \\
awk '{{OFS="\\t"; print ">"$1"\\n"$10}}' - \\
  > {fasta_file}""".format(
                    input=input,
                    fasta_file=fasta_file
                )
                current_jobs.append(
                    Job(
                        [input],
                        [],
                        [["blast", "module_samtools"]],
                        command=command
                    )
                )

                # run blast
                command = """\
blastn \\
  -query {fasta_file} \\
  -db nt \\
  -out {result_file} \\
  -perc_identity 80 \\
  -num_descriptions 1 \\
  -num_alignments 1""".format(
                    fasta_file=fasta_file,
                    result_file=result_file
                )
                current_jobs.append(
                    Job(
                        [],
                        [],
                        [["blast", "module_blast"]],
                        command=command
                    )
                )

                # filter and format the result to only have the sorted number of match and the species
                command = """\
grep ">" {result_file} | \\
awk ' {{ print $2 "_" $3}} ' | \\
sort | \\
uniq -c | \\
sort -n -r | \\
head -20 > {output} && true""".format(
                    result_file=result_file,
                    output=output
                )
                current_jobs.append(
                    Job(
                        [],
                        [output],
                        [],
                        command=command
                    )
                )
            else:
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
                    nb_blast_to_do=nb_blast_to_do)
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

                command = """
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

                command = """
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
                name="blast." + readset.name + ".blast." + self.run_id + "." + str(self.lane_number)
                samples = [readset.sample]
            )
            jobs.append(job)

        self.add_copy_job_inputs(jobs)
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

        for readset in self.readsets:
            output_dir = os.path.join(os.path.dirname(readset.fastq1), "qc")
            region_name = readset.name + "_" + readset.sample_number + "_L00" + readset.lane

            file1 = readset.fastq1
            file2 = readset.fastq2
            type = "FASTQ"
            if readset.bam:
                file1 = readset.bam + ".bam"
                file2 = None
                type = "BAM"

            jobs.append(
                concat_jobs([
                    bash.mkdir(output_dir),
                    bvatools.readsqc(
                        file1,
                        file2,
                        type,
                        region_name,
                        output_dir
                )],
                name="qc." + readset.name + ".qc." + self.run_id + "." + str(self.lane_number),
                samples=[readset.sample]
            ))

        self.add_copy_job_inputs(jobs)
        return self.throttle_jobs(jobs)

    def fastqc(self):
        """
        """
        jobs = []
        for sample in self.samples:
            fastqc_directory = os.path.join(self.output_dir, "Unaligned." + readset.lane, 'Project_' + readset.project, 'Sample_' + readset.name, "fastqc.R1")  
            input = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.bam")
            output = os.path.join(fastqc_directory, sample.name + ".sorted.dup_fastqc.zip")

            jobs.append(
                concat_jobs([
                    bash.mkdir(
                        fastqc_directory,
                        remove=True
                    ),
                    adapter_job,
                    fastqc.fastqc(
                        input,
                        None,
                        output,
                        adapter_file=None,
                        use_tmp=True
                )],
                name="fastqc." + sample.name,
                samples=[sample]
            ))

    def md5(self):
        """
        Create md5 checksum files for the fastq, bam and bai using the system 'md5sum'
        util.

        One checksum file is created for each file.
        """
        jobs = []
        for readset in self.readsets:
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

            jobs.append(
                concat_jobs(
                    current_jobs
                ),
                name="md5." + readset.name + ".md5." + self.run_id + "." + str(self.lane_number),
                samples=[readset.sample]
            )

        if config.param('md5', 'one_job', required=False, type="boolean"):
            job = concat_jobs(
                jobs,
                name="md5." + self.run_id + "." + str(self.lane_number),
                samples=self.readsets
            )
            self.add_copy_job_inputs([job])
            return [job]
        else:
            self.add_copy_job_inputs(jobs)
            return jobs

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

        # Notification
        output1 = os.path.join(self.output_dir, "notificationProcessingComplete." + str(self.lane_number) + ".out")
        output2 = os.path.join(self.output_dir, "notificationCopyStart." + str(self.lane_number) + ".out")

        notification_command = config.param('copy', 'notification_command', required=False)
        if notification_command:
            job = Job(
                inputs,
                [output1, output2],
                name="start_copy_notification." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            )
            job.command = notification_command.format(
                technology=config.param('copy', 'technology'),
                output_dir=self.output_dir,
                run_id=self.run_id,
                output1=output1,
                output2=output2,
                lane_number=self.lane_number
            )
            jobs_to_concat.append(job)

        # Actual copy
        full_destination_folder = config.param('copy', 'destination_folder', type="dirpath") + os.path.basename(self.run_dir)
        output = os.path.join(full_destination_folder, "copyCompleted." + str(self.lane_number) + ".out")

        exclude_bam = config.param('copy', 'exclude_bam', required=False, type='boolean')
        exclude_fastq_with_bam = config.param('copy', 'exclude_fastq_with_bam', required=False, type='boolean')
        if exclude_bam and exclude_fastq_with_bam:
            log.warn("Excluding both BAM and fastq files")

        excluded_files = []

        if exclude_bam or exclude_fastq_with_bam:
            for readset in [readset for readset in self.readsets if readset.bam]:
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
                lane_number=self.lane_number,
                run_id=self.run_id,
                source=self.run_dir,
                run_name=os.path.basename(self.run_dir)
            )
            jobs_to_concat.append(
                Job(
                    inputs,
                    [output],
                    command=copy_command_run_folder,
                    samples=self.samples
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
                samples=self.samples
            )
        )
        jobs_to_concat.append(
            Job(
                command="touch " + output,
                samples=self.samples
            )
        )

        job = concat_jobs(
            jobs_to_concat,
            name="copy." + self.run_id + "." + str(self.lane_number)
        )

        return [job]

    def end_copy_notification(self):
        """
        Send an optional notification to notify that the copy is finished.

        The command used is in the configuration file. This step is skipped when no
        command is provided.
        """
        jobs = []

        full_destination_folder = config.param('copy', 'destination_folder', type="dirpath") + os.path.basename(
            self.run_dir)
        input = os.path.join(full_destination_folder, "copyCompleted." + str(self.lane_number) + ".out")
        output = os.path.join(full_destination_folder, "notificationAssociation." + str(self.lane_number) + ".out")

        notification_command = config.param('end_copy_notification', 'notification_command', required=False)
        if notification_command:
            job = Job(
                [input],
                [output],
                name="end_copy_notification." + self.run_id + "." + str(self.lane_number)
            )
            job.command = notification_command.format(
                technology=config.param('end_copy_notification', 'technology'),
                output_dir=self.output_dir,
                run_name=os.path.basename(self.run_dir),
                run_id=self.run_id,
                output=output,
                lane_number=self.lane_number
            )
            job.samples = self.samples
            jobs.append(job)

        return jobs

    #
    # Utility methods
    #

    def add_copy_job_inputs(self, jobs):
        for job in jobs:
            # we first remove dependencies of the current job, since we will have a dependency on that job
            self.copy_job_inputs = [item for item in self.copy_job_inputs if item not in job.input_files]
            self.copy_job_inputs.extend(job.output_files)

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
        for read in self.read_infos if read.is_index:
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
        instrument_file = config.param('DEFAULT', 'instrument_list_file', type='filepath', required='false')
        if not (instrument_file and os.path.isfile(instrument_file)):
            instrument_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'instrument_list.csv')

        return subprocess.check_output("grep -m1 '"+instrument+"' %s | awk -F',' '{print $3}'" % instrument_file, shell=True).strip()

    def get_instrument(self):
        """
        Parse the RunInfo.xml file for the instrument name the run has been running on
        """
        return Xml.parse(os.path.join(self.run_dir, "RunInfo.xml")).getroot().find('Run').find('Instrument').text

    def validate_barcodes(self):
        """
        Validate all index sequences against each other to ensure they aren't in collision according to the chosen
        number of mismatches parameter.
        """
        min_allowed_distance = (2 * self.number_of_mismatches) + 1

        validated_indexes = []
        collisions = []

        for readset in self.readsets:
            current_index = readset.index.replace('-', '')

            for candidate_index in validated_indexes:
                if distance(current_index, candidate_index) < min_allowed_distance:
                    collisions.append("'" + current_index + "' and '" + candidate_index + "'")
            validated_indexes.append(current_index)

        if len(collisions) > 0:
            _raise(SanitycheckError("Barcode collisions: " + ";".join(collisions)))

    def get_mask(self):
        """
        Returns a BCL2FASTQ friendly mask of the reads cycles.

        The mask is calculated using:
            - first base and last base of index;
            - the index length in the sample sheet;
            - the number of index cycles on the sequencer;
        """
        mask = ""
        index_lengths = self.get_smallest_index_length()
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

    def get_indexcycles(self):
        for read in self.read_infos if read.is_index:
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
        instrument_file = config.param('DEFAULT', 'instrument_list_file', type='filepath', required='false')
        if not (instrument_file and os.path.isfile(instrument_file)):
            instrument_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'instrument_list.csv')

        return subprocess.check_output("grep -m1 '"+instrument+"' %s | awk -F',' '{print $3}'" % instrument_file, shell=True).strip()

    def get_instrument(self):
        """
        Parse the RunInfo.xml file for the instrument name the run has been running on
        """
        return Xml.parse(os.path.join(self.run_dir, "RunInfo.xml")).getroot().find('Run').find('Instrument').text

    def generate_illumina_lane_sample_sheet(self):
        """
        Create a sample sheet to use with the BCL2FASTQ software.

        Only the samples of the chosen lane will be in the file.
        The sample indexes are trimmed according to the mask used.
        """

        if self.protocol == "nanuq":
            jobs = self.generate_nanuq_sample_sheet()
        else:
            jobs = self.generate_clarity_sample_sheet()
        
        return jobs

    def generate_nanuq_sample_sheet(self):
        """
        Create a sample sheet for bcl2fastq from the Nanuq data
        """

        read_masks = self.mask.split(",")
        has_single_index = self.has_single_index()

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
            config.param('DEFAULT', 'casava_sample_sheet_prefix') + str(self.lane_number) + ".csv"
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

        for readset in self.readsets:
            index_to_use = ""

            if len(readset.index) > 0 and len(self.readsets) > 1:
                indexes = readset.index.split("-")
                nb_index = len(indexes)

                if has_single_index:
                    # we have a mixed of index in the sample, there are samples with 1 or 2 index,
                    # ignore the second index in the samplesheet
                    nb_index = 1

                for i in range(0, nb_index):
                    nb_ignored_leading_bases = 0
                    nb_of_index_bases = 0

                    m = re.match("(n\d+)?(I\d+)(n\d+)?", read_masks[i + 1])
                    if m:
                        if m.group(1):
                            nb_ignored_leading_bases = int(m.group(1)[1:])
                        if m.group(2):
                            nb_of_index_bases = int(m.group(2)[1:])

                    # remove ignored leading bases and trim index to smallest lane index
                    index = indexes[i][nb_ignored_leading_bases:nb_ignored_leading_bases + nb_of_index_bases]

                    if i > 0 and len(index) > 0:
                        index_to_use += "-"
                    index_to_use += index

            readset._index = index_to_use if len(index_to_use) > 0 else "NoIndex"
            index_array = index_to_use.split("-")

            fastq_file_pattern = os.path.join(
                self.output_dir,
                "Unaligned." + readset.lane,
                'Project_' + readset.project,
                'Sample_' + readset.name,
                readset.name + '_S' + readset.sample_number + '_L00' + readset.lane +
                '_R{read_number}_001.fastq.gz'
            )
            readset.fastq1 = fastq_file_pattern.format(read_number=1)
            readset.fastq2 = fastq_file_pattern.format(read_number=2) if readset.run_type == "PAIRED_END" else None
            readset.index_fastq1 = re.sub("_R1_", "_I1_", readset.fastq1)
            readset.index_fastq2 = re.sub("_R2_", "_I2_", readset.fastq2) if readset.run_type == "PAIRED_END" and readset.index_type == "DUALINDEX" else None

            csv_dict = {
                "FCID": readset.flow_cell,
                "Lane": self.lane_number,
                "Sample_ID": "Sample_" + readset.name,
                "Sample_Name": readset.name,
                "SampleRef": "",
                "Index": index_array[0],
                "Index2": index_array[1] if len(index_array) > 1 else "",
                "Description": readset.description,
                "Control": readset.control,
                "Recipe": readset.recipe,
                "Operator": readset.operator,
                "Sample_Project": "Project_" + readset.project
            }
            writer.writerow(csv_dict)

    def generate_clarity_sample_sheet(self):
        """
        Create a sample sheet for bcl2fastq from the Clarity data
        """

        mkdir_job = bash.mkdir(os.path.join(self.output_dir, "Unaligned." + readset.lane, 'Project_' + readset.project, 'Sample_' + readset.name))

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
            config.param('DEFAULT', 'casava_sample_sheet_prefix') + str(self.lane_number) + ".csv"
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

        mask = self.mask
        index1cycles = self.index1cycles
        index2cycles = self.index2cycles

        overmask = ""
        overindex1 = None
        overindex2 = None

        dualindex_demultiplexing = False    # This is not the sequencing demultiplexing flag, but really the index demultiplexing flag
        merge_undetermined = False
        generate_umi = False

        # If only one library on the lane
        if len(self.readsets) == 1:
            merge_undetermined = True
        merge = config.param('fastq', 'merge_undetermined') if config.param('fastq', 'merge_undetermined', required=False, type='boolean') else merge_undetermined

        # barcode validation
        if re.search("I", mask):
            self.validate_barcodes()

        # IDT - UMI9 in index2
        if re.search(",I17", mask):
            generate_umi = True
            overmask = re.sub(",I17,", ",I8n*,", mask)
            overindex1=8
            overindex2=8

        # HaloPlex - UMI8 in index2 alone
        bcl2fastq_extra_option = ""
        if sum(1 for readset in [readset for readset in self.readsets if re.search("HaloPlex", readset.index)]) > 0:
            if "DUALINDEX" in set([readset.index_type for readset in self.readsets]) :
                _raise(SanityCheckError("HaloPlex libraries cannot be mixed with DUAL INDEX libraries"))
            overmask=re.sub(",I8,I10,", ",I8,Y10,", mask)
            overindex1=8
            overindex2=0
            bcl2fastq_extra_option="--mask-short-adapter-reads 10"

        # If SINGLEINDEX only
        if "DUALINDEX" not in set([readset.index_type for readset in self.readsets]):
            if re.search("I", mask.split(",")[2]):
                split_mask = mask.split(",") if overmask == "" else overmask.split(",")
                overmask=','.join(split_mask[0], split_mask[1], "n*", split_mask[3])
                overindex1=8
                overindex2=0

        final_mask = mask if overmask == "" else overmask
        final_index1 = self.index1cycles if overindex1 is None else overindex1
        final_index2 = self.index2cycles if overindex2 is None else overindex2

        mask = config.param('fastq', 'overmask') if config.param('fastq', 'overmask', required=False, type='string') else final_mask
        self._index1cycles = config.param('fastq', 'overindex1') if config.param('fastq', 'overindex1', required=False, type='int') else final_index1
        self._index2cycles = config.param('fastq', 'overindex2') if config.param('fastq', 'overindex2', required=False, type='int') else final_index2

        # If the second index exists
        if self.index2cycles != 0:
            dualindex_demultiplexing = True

        # In case of HaloPlex-like masks, R2 is actually I2 (while R3 is R2),
        # Proceed to dualindex demultiplexing
        if ''.join(i for i in mask if not i.isdigit()) == "Y,I,Y,Y":
            dualindex_demultiplexing = True

        output_dir = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number))
        casava_sample_sheet = os.path.join(self.output_dir, "casavasheet." + str(self.lane_number) + ".indexed.csv")

        count = 0
        index_per_sample = []
        for readset in self.readsets:
            count += 1
            readset_indexes = self.get_index(readset)

            indexes = []
            for readset_index in readset_indexes:
                    indexes.append({
                        readset.description : [index['INDEX1'])

                index_array = [readset_index['INDEX1'], readset_index['INDEX2']]
                indexes.append({readset.description : index_array})

                fastq_file_pattern = os.path.join(
                    self.output_dir,
                    "Unaligned." + readset.lane,
                    'Project_' + readset.project,
                    'Sample_' + readset.name,
                    readset.name + '_S' + readset.sample_number + '_L00' + readset.lane +
                    '_R{read_number}_001.fastq.gz'
                )
                readset.fastq1 = fastq_file_pattern.format(read_number=1)
                readset.fastq2 = fastq_file_pattern.format(read_number=2) if readset.run_type == "PAIRED_END" else None
                readset.index_fastq1 = re.sub("_R1_", "_I1_", readset.fastq1)
                readset.index_fastq2 = re.sub("_R2_", "_I2_", readset.fastq2) if readset.run_type == "PAIRED_END" and readset.index_type == "DUALINDEX" else None

                csv_dict = {
                    "FCID": readset.flow_cell,
                    "Lane": self.lane_number,
                    "Sample_ID": "Sample_" + readset_index['BC2LFASTQ_NAME'],
                    "Sample_Name": readset_index['BC2LFASTQ_NAME'],
                    "SampleRef": "",
                    "Index": index_array[0],
                    "Index2": index_array[1],
                    "Description": readset.description + ' - ' + readset.lib_type + ' - ' + readset.library_source
                    "Control": readset.control,
                    "Recipe": readset.recipe,
                    "Operator": readset.operator,
                    "Sample_Project": "Project_" + readset.project
                }
                writer.writerow(csv_dict)
            index_per_sample.append({readset.name : indexes})
        self._indexes_from_lims = index_per_sample

    def generate_fastq_outputs(self, merge=False):
        bcl2fastq_outputs = []
        count = 0

        merge_undetermined = False
        generate_umi = False

        # If only one library on the lane
        if len(self.readsets) == 1:
            merge_undetermined = True
        merge = config.param('fastq', 'merge_undetermined') if config.param('fastq', 'merge_undetermined', required=False, type='boolean') else merge_undetermined

        output_dir = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number))
        casava_sample_sheet = os.path.join(self.output_dir, "casavasheet." + str(self.lane_number) + ".indexed.csv")

        # IDT - UMI9 in index2
        if re.search(",I17", mask):
            generate_umi = True
            overmask = re.sub(",I17,", ",I8n*,", mask)
            overindex1=8
            overindex2=8

        if generate_umi:
            output_dir_noindex = output_dir + ".noindex")
            casava_sample_sheet_noindex = re.sub(".indexed.", ".noindex.", casava_sample_sheet) 

        for readset in self.readsets:
            count += 1

            # If 10X libraries : 4 indexes per sample
            if re.search("tenX", readset.library_type):
                fastq1 = os.path.basename(readset.fastq1)
                counta = count
                bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_A_" + counta + "_", fastq1))
                countb = count + 1
                bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_B_" + countb + "_", fastq1))
                countc = count + 2
                bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_C_" + countc + "_", fastq1))
                countd = count + 3
                bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_D_" + countd + "_", fastq1))
                # If True, then merge the 'Undetermined' reads
                if merge:
                    bcl2fastq_outputs.append(re.sub(readset.name, "Undetermined_S0", fastq1))

                # Add the fastq of first index
                idx_fastq1 = re.sub("R1", "I1", fastq1)
                bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_A_" + counta + "_", idx_fastq1))
                bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_B_" + countb + "_", idx_fastq1))
                bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_C_" + countc + "_", idx_fastq1))
                bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_D_" + countd + "_", idx_fastq1))
                # If True, then merge the 'Undetermined' reads
                if merge:
                    bcl2fastq_outputs.append(re.sub(readset.name, "Undetermined_S0", idx_fastq1))

                # For paired-end sequencing, do not forget the fastq of the reverse reads
                if readset.run_type == "PAIRED_END" :
                    fastq2 = os.path.basename(readset.fastq2)
                    bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_A_" + counta + "_", fastq2))
                    bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_B_" + countb + "_", fastq2))
                    bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_C_" + countc + "_", fastq2))
                    bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_D_" + countd + "_", fastq2))
                    # If True, then merge the 'Undetermined' reads
                    if merge:
                        bcl2fastq_outputs.append(re.sub(readset.name, "Undetermined_S0", fastq2))

                    # For dual index demultiplexing, do not forget the fastq of the second index
                    idx_fastq2 = re.sub("R2", "I2", fastq2)
                    bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_A_" + counta + "_", idx_fastq2))
                    bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_B_" + countb + "_", idx_fastq2))
                    bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_C_" + countc + "_", idx_fastq2))
                    bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_D_" + countd + "_", idx_fastq2))
                    # If True, then merge the 'Undetermined' reads
                    if merge:
                        bcl2fastq_outputs.append(re.sub(readset.name, "Undetermined_S0", idx_fastq2))

                count = countd

            # not a 10X library : 1 index per sample
            else:
                fastq1 = os.path.basename(readset.fastq1)
                bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_" + count + "_", fastq1))
                # If ask to merge the Undeternined reads
                if merge:
                    bcl2fastq_outputs.append(re.sub(readset.name, "Undetermined_S0", fastq1))

                idx_fastq1 = re.sub("R1", "I1", fastq1)
                bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_" + count + "_", idx_fastq1))
                # If ask to merge the Undeternined reads
                if merge:
                    bcl2fastq_outputs.append(re.sub(readset.name, "Undetermined_S0", idx_fastq1))

                # For paired-end sequencing, do not forget the fastq of the reverse reads
                if readset.run_type == "PAIRED_END" :
                    fastq2 = os.path.basename(readset.fastq2)
                    bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_" + count + "_", fastq2))
                    # If True, then merge the 'Undetermined' reads
                    if merge:
                        bcl2fastq_outputs.append(re.sub(readset.name, "Undetermined_S0", fastq2))

                    # For dual index multiplexing, do not forget the fastq of the second index
                    idx_fastq2 = re.sub("R2", "I2", fastq2)
                    bcl2fastq_outputs.append(re.sub(readset.name, readset.name + "_" + count + "_", idx_fastq2))
                    # If True, then merge the 'Undetermined' reads
                    if merge:
                        bcl2fastq_outputs.append(re.sub(readset.name, "Undetermined_S0", idx_fastq2))

        return bcl2fastq_outputs

    def get_index(self, readset):
        """
        """

        indexes = []
        index1 = ""
        index2 = ""
        index1seq = ""
        index2seq = ""

        index_file = config.param('DEFAULT', 'index_settings_file', type='filepath', required='false')
        if not (index_file and os.path.isfile(index_file)):
            index_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_settings_format.txt')

        index_seq_pattern = "grep '%s,' %s | head -n1 | awk -F',' '{print $%d}'"

        if readset.library_type == "tenX_sc_RNA_v2" or readset.library_type == "tenX_DNA_v2" :
            index1 = readset.description

            for idx, char in enumerate(['_A', '_B', '_C', '_D']):
                index1seq = subprocess.check_output(index_seq_pattern % (index1, index_file, idx+2), shell=True).strip()
                [actual_index1seq, actual_index2seq, adapteri7, adapteri5] = self.sub_get_index(readset, index1seq, index2seq)
                indexes.append({
                    'BC2LFASTQ_NAME': readset.name + char,
                    'INDEX1': actual_index1seq,
                    'INDEX2': actual_index2seq,
                    'ADAPTERi7' : adapteri7,
                    'ADAPTERi5' : adapteri5
                })

        elif readset.library_type == "tenX_sc_RNA_v1" :
            index2 = readset.description 

            for idx, char in enumerate(['_A', '_B', '_C', '_D']):
                index2seq = subprocess.check_output(index_seq_pattern % (index2, index_file, idx+2), shell=True).strip()
                [actual_index1seq, actual_index2seq, adapteri7, adapteri5] = self.sub_get_index(readset, index1seq, index2seq)
                indexes.append({
                    'BC2LFASTQ_NAME': readset.name + char,
                    'INDEX1': actual_index1seq,
                    'INDEX2': actual_index2seq,
                    'ADAPTERi7' : adapteri7,
                    'ADAPTERi5' : adapteri5
                })

        elif re.search("-", readset.description):
            index1 = readset.description.split("-")[0]
            index2 = readset.description.split("-")[1]

            index1seq = subprocess.check_output(index_seq_pattern % (index1, index_file, 2), shell=True).strip()
            index2seq = subprocess.check_output(index_seq_pattern % (index2, index_file, 2), shell=True).strip()

            [actual_index1seq, actual_index2seq, adapteri7, adapteri5] = self.sub_get_index(readset, index1seq, index2seq)
            indexes.append({
                'BC2LFASTQ_NAME': readset.name,
                'INDEX1': actual_index1seq,
                'INDEX2': actual_index2seq,
                'ADAPTERi7' : adapteri7,
                'ADAPTERi5' : adapteri5
            })

        else:
            index1 = readset.description
            
            index1seq = subprocess.check_output(index_seq_pattern % (index1, index_file, 2), shell=True).strip()

            [actual_index1seq, actual_index2seq, adapteri7, adapteri5] = self.sub_get_index(readset, index1seq, index2seq)
            indexes.append({
                'BC2LFASTQ_NAME': readset.name,
                'INDEX1': actual_index1seq,
                'INDEX2': actual_index2seq,
                'ADAPTERi7' : adapteri7,
                'ADAPTERi5' : adapteri5
            })

        return indexes

    def sub_get_index(self, readset, index1seq, index2seq)
        """
        """

        index_file = config.param('DEFAULT', 'index_settings_file', type='filepath', required='false')
        if not (index_file and os.path.isfile(index_file)):
            index_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_settings_format.txt')

        actual_index1seq='';
        actual_index2seq='';

        seqtype = self.seqtype

        primer_seq_pattern = "grep -A8 '%s' %s | grep '_IDX_' | awk -F':' '{print $2}' | tr -d \"35\'\- \" | awk -F',' '{print $%d}'"

        index1_primer = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index 1") % (seqtype, index_file, 1), shell=True).strip()
        index1_primeroffset = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index 1") % (seqtype, index_file, 2), shell=True).strip()
        index2_primer = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index 2") % (seqtype, index_file, 1), shell=True).strip()
        index2_primeroffset = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index 2") % (seqtype, index_file, 2), shell=True).strip()

        main_seq_pattern = "grep -A4 '^%s:' %s | grep -E \"^3'|^5'\" | head -n 1 | sed \"s/5'//g\"  | sed \"s/3'//g\" | tr -d \" '\-\" | grep %s"
        actual_seq_pattern = "echo %s | awk -F\"%s\" '{print $%d}' | sed \"%s\" | cut -c %d"

        present = subprocess.check_output(main_seq_pattern.replace("grep %s", "grep -c %s") % (readset.library_type, index_file, index1_primer), shell=True).strip()

        if present == "1" :

            main_seq = subprocess.check_output(main_seq_pattern % (readset.library_type, index_file, index1_primer), shell=True).strip()

            actual_index1seq = subprocess.check_output(actual_seq_pattern % (main_seq, index1_primer, 2, "s/\[i7\]/index1seq/g", int(index1_primeroffset) + 1 - int(index1_primeroffset) + int(self.index1cycles)), shell=True).strip()

            if readset.index_type == "DUALINDEX" :
                main_seq = subprocess.check_output(main_seq_pattern.replace("| head -n 1 |", "|") % (readset.library_type, index_file, index2_primer), shell=True).strip()

                if seqtype == "hiseqx" or seqtype == "hiseq4000" or seqtype == "iSeq" :
                    actual_index2seq = subprocess.check_output(actual_seq_pattern.replace("| cut -c", "| rev | cut -c") % (main_seq, index2_primer, 1, "s/\[i5c\]/$(echo "+index2seq+" | tr 'ATGC' 'TACG' )/g", int(index2_primeroffset) + 1 - int($index2_primeroffset) + int(self.index2cycles)), shell=True).strip()
                else :
                    actual_index2seq = subprocess.check_output(actual_seq_pattern % (main_seq, index2_primer, 2, "s/\[i5\]/index2seq/g", int(index2_primeroffset) + 1 - int(index2_primeroffset) + int(self.index2cycles)), shell=True).strip()

        else :
            indexn1_primer = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index N1").replace("_DIGIT_", "1"), shell=True).strip()
            indexn1_primeroffset = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index N1").replace("_DIGIT_", "2"), shell=True).strip()
            indexn2_primer = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index N2").replace("_DIGIT_", "1"), shell=True).strip()
            indexn2_primeroffset = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index N2").replace("_DIGIT_", "2"), shell=True).strip()

            main_seq = subprocess.check_output(main_seq_pattern % (readset.library_type, index_file, indexn1_primer), shell=True).strip()

            actual_index1seq = subprocess.check_output(actual_seq_pattern % (main_seq, indexn1_primer, 2, "s/\[i7\]/index1seq/g", int($indexn1_primeroffset) + 1 - int(indexn1_primeroffset) + int(self.index1cycles)), shell=True).strip()

            if readset.index_type == "DUALINDEX" :
                main_seq = subprocess.check_output(main_seq_pattern.replace("| head -n 1 |", "|") % (readset.library_type, index_file, indexn2_primer), shell=True).strip()

                if seqtype == "hiseqx" or seqtype == "hiseq4000" or seqtype == "iSeq" :
                    actual_index2seq = subprocess.check_output(actual_seq_pattern.replace("| cut -c", "| rev | cut -c") % (main_seq, indexn2_primer, 1, "s/\[i5c\]/$(echo "+index2seq+" | tr 'ATGC' 'TACG' )/g", int(indexn2_primeroffset) + 1 - int(indexn2_primeroffset) + int(self.self.index2cycles)), shell=True).strip()
                else :
                    actual_index2seq = subprocess.check_output(actual_seq_pattern % (main_seq, indexn2_primer, 2, "s/\[i5\]/index2seq/g", int(indexn2_primeroffset) + 1 - int(indexn2_primeroffset) + int(self.index2cycles)), shell=True).strip()

        main_seq = subprocess.check_output(main_seq_pattern % (readset.library_type, index_file, "\[i7\]"), shell=True).strip()
        adapteri7=subprocess.check_output("echo \"%s\" | awk -F \"\\\[i7\\\]\" '{print $1}' | awk -F \"\\\]\" '{print $NF}'" % (main_seq), shell=True).strip()

        main_seq = subprocess.check_output(main_seq_pattern.replace("| head -n 1 |", "|") % (readset.library_type, index_file, "\[i5c\]"), shell=True).strip()
        adapteri5=subprocess.check_output("echo \"%s\" | awk -F \"\\\[i5c\\\]\" '{print $2}' | awk -F \"\\\[\" '{print $1}' | rev" % (main_seq), shell=True).strip()

        return [actual_index1seq, actual_index2seq, adapteri7, adapteri5]

    def has_single_index(self):
        """ 
        Returns True when there is at least one sample on the lane that doesn't use double-indexing or we only have
        one read of indexes.
        """
        return len([readset for readset in self.readsets if ("-" not in readset.index)]) > 0 or\
               len([read for read in self.read_infos if read.is_index]) < 2

    def get_smallest_index_length(self):
        """
        Returns a list (for each index read of the run) of the minimum between the number of index cycle on the
        sequencer and all the index lengths.
        """
        run_index_lengths = [r.nb_cycles for r in self.read_infos if r.is_index] # from RunInfo

        if len(run_index_lengths) == 0 and len(self.readsets) > 1:
            _raise(SanitycheckError("Multiple samples on a lane, but no indexes were read from the sequencer."))

        # loop on all index reads, to compare with samples index length
        for i in range(0, len(run_index_lengths)):
            min_sample_index_length = 0
            try:
                min_sample_index_length = min(len(readset.index.split("-")[i])
                                              for readset in self.readsets
                                              if (len(readset.index.split("-")) > i and len(readset.index.split("-")[i]) > 0))
            except ValueError:
                pass  # we don't have a sample with this Ith index read, use the 0 already set

            empty_index_list = [ readset for readset in self.readsets if (len(readset.index.split("-")) <= i or len(readset.index.split("-")[i]) == 0) ]
        
            if len(empty_index_list):
                # we have samples without this Ith index read, so we skip it
                min_sample_index_length = 0

            run_index_lengths[i] = min(min_sample_index_length, run_index_lengths[i])

        return run_index_lengths

    def interop_getrawreads(self)
        """
        Parse the InterOp .bin files (within the run folder) to retrieve the # raw reads
        """
        run_metrics = py_interop_run_metrics.run_metrics()

        valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
        py_interop_run_metrics.list_index_metrics_to_load(valid_to_load)
        run_metrics.read(self.run_dir, valid_to_load)

        summary = py_interop_summary.index_flowcell_summary()
        py_interop_summary.summarize_index_metrics(run_metrics, summary)

        lane_summary = summary.at(self.lane_number - 1)
        
        return lane_summary.total_reads()

    def parse_run_info_file(self):
        """
        Parse the RunInfo.xml file of the run and returns the list of RunInfoRead objects
        """
        reads = Xml.parse(os.path.join(self.run_dir, "RunInfo.xml")).getroot().find('Run').find('Reads')
        return [ RunInfoRead(int(r.get("Number")), int(r.get("NumCycles")), r.get("IsIndexedRead") == "Y") for r in reads.iter('Read') ]

    def load_readsets(self):
        """
        Download the sample sheets if required or asked for; call the load of these files and return a list of
        readsets.
        """

        if self.protocol == "nanuq":
            # Casava sheet download
            if not self.args.casava_sheet_file or self.args.force_download:
                if not os.path.exists(self.casava_sheet_file) or self.args.force_download:
                    command = config.param('DEFAULT', 'fetch_casava_sheet_command').format(
                        output_directory=self.output_dir,
                        run_id=self.run_id,
                        filename=self.casava_sheet_file
                    )
                    log.info(command)
                    return_code = subprocess.call(
                        command,
                        shell=True
                    )
                    if return_code != 0:
                        _raise(SanitycheckError("Unable to download the Casava Sheet."))

            # Nanuq readset file download
            if not self.args.readsets or self.args.force_download:
                if not os.path.exists(self.nanuq_readset_file) or self.args.force_download:
                    command = config.param('DEFAULT', 'fetch_nanuq_sheet_command').format(
                        output_directory=self.output_dir,
                        run_id=self.run_id,
                        filename=self.nanuq_readset_file
                    )
                    return_code = subprocess.call(
                        command,
                        shell=True
                    )
                    if return_code != 0:
                        _raise(SanitycheckError("Unable to download the Nanuq readset file."))

            return parse_illumina_raw_readset_files(
                self.output_dir,
                self.run_dir,
                "PAIRED_END" if self.is_paired_end else "SINGLE_END",
                self.nanuq_readset_file,
                self.casava_sheet_file,
                self.args.lane_number,
                config.param('DEFAULT', 'genomes_home', type="dirpath"),
                self.get_sequencer_minimum_read_length(),
                self.protocol
            )
        
        else:
            return parse_illumina_raw_readset_files(
                self.output_dir,
                self.run_dir,
                "PAIRED_END" if self.is_paired_end else "SINGLE_END",
                self.clarity_event_file,
                None,
                self.args.lane_number,
                config.param('DEFAULT', 'genomes_home', type="dirpath"),
                self.get_sequencer_minimum_read_length(),
                self.protocol
            )

    def submit_jobs(self):
        super(IlluminaRunProcessing, self).submit_jobs()

    def throttle_jobs(self, jobs):
        """
        Group jobs of the same task (same name prefix) if they exceed the configured threshold number.
        """

        max_jobs_per_step = config.param('default', 'max_jobs_per_step', required=False, type="int")
        jobs_by_name = collections.OrderedDict()
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
                        merged_jobs.append(concat_jobs(
                            current_jobs[x * number_task_by_job:min((x + 1) * number_task_by_job, len(current_jobs))],
                            job_name + "." + str(x + 1) + "." + self.run_id + "." + str(self.lane_number)))
                reply.extend(merged_jobs)
            else:
                reply.extend(current_jobs)
        return reply

    @property
    def steps(self):
        return [
            [
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
                self.copy,
                self.end_copy_notification
            ],[
                self.index,
                self.fastq,
                self.align,
                self.picard_mark_duplicates,
                self.metrics,
                self.blast,
                self.qc_graphs,
                self.md5,
                self.copy,
                self.end_copy_notification
            ]
        ]

def distance(str1, str2):
    """
    Returns the hamming distance. http://code.activestate.com/recipes/499304-hamming-distance/#c2
    """
    return sum(itertools.imap(unicode.__ne__, str1, str2))

if __name__ == '__main__':
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        pipeline = IlluminaRunProcessing(protocol=['clarity', 'nanuq'])
