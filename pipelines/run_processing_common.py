#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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
import argparse
import logging
import os
import re
import socket
import string
import math
import csv
import sys
import subprocess
import json
from collections import OrderedDict

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs 

from bfx import bvatools
from bfx import picard
from bfx import fastp
from bfx import fastqc
from bfx import tools
from bfx import bash_cmd as bash

from pipelines import common

log = logging.getLogger(__name__)

class RunProcessing(common.MUGQICPipeline):

    def __init__(self, protocol=None):
        self._protocol=protocol
        self.copy_job_inputs = {}
        self.argparser.add_argument("-r", "--readsets", help="Sample sheet for the MGI run to process (mandatory)", type=argparse.FileType('r'), required=False)
        self.argparser.add_argument("-d", "--run", help="Run directory (mandatory)", required=False, dest="run_dir")
        self.argparser.add_argument("--run-id", help="Run ID. Default is parsed from the run folder", required=False, dest="run_id")
        self.argparser.add_argument("--flowcell-id", help="Flowcell ID. Default is parsed from the run folder", required=False, dest="flowcell_id")
        self.argparser.add_argument("--lane", help="Lane number (to only process the given lane)", type=int, required=False, dest="lane_number")
        self.argparser.add_argument("-x", help="First index base to use for demultiplexing (inclusive). The index from the sample sheet will be adjusted according to that value.", type=int, required=False, dest="first_index")
        self.argparser.add_argument("-y", help="Last index base to use for demultiplexing (inclusive)", type=int, required=False, dest="last_index")
        self.argparser.add_argument("-m", help="Number of index mistmaches allowed for demultiplexing (default 1). Barcode collisions are always checked.", type=int, required=False, dest="number_of_mismatches")
        self.argparser.add_argument("--allow-barcode-collision", help="Allow barcode collision by not comparing barcode sequences to each other (usually decreases the demultiplexing efficiency).", action="store_true", required=False, dest="allow_barcode_collision")

        super(RunProcessing, self).__init__(protocol)

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = {}
            if not hasattr(self, "_mask"):
                self._mask = {} 
            for lane in self.lanes:
                self._readsets[lane] = self.load_readsets(lane)
                if not int(self.index1cycles[lane]) + int(self.index2cycles[lane]) == 0:
                    self._mask[lane] = self.get_mask(lane)
                    self.generate_lane_sample_sheet(lane)
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
                self._lanes = [str(self.lane_number)]
            else:
                self._lanes = [lane for lane in list(set([line['Position'].split(":")[0] for line in csv.DictReader(open(self.readset_file, 'r'), delimiter='\t', quotechar='"')]))]
            for lane in self._lanes:
                self.copy_job_inputs[lane] = []
        return self._lanes

    @property
    def is_paired_end(self):
        if not hasattr(self, "_is_paired_end"):
            self._is_paired_end = {} 
            for lane in self.lanes:
                if self.read2cycles[lane]:
                    self._is_paired_end[lane] = True 
                else:
                    self._is_paired_end[lane] = False
        return self._is_paired_end

    @property
    def read1cycles(self):
        if not hasattr(self, "_read1cycles"):
            self._read1cycles = {} 
            for lane in self.lanes:
                self._read1cycles[lane] = self.get_read1cycles(lane)
        return self._read1cycles

    @property
    def read2cycles(self):
        if not hasattr(self, "_read2cycles"):
            self._read2cycles = {} 
            for lane in self.lanes:
                self._read2cycles[lane] = self.get_read2cycles(lane)
        return self._read2cycles

    @property
    def is_dual_index(self):
        if not hasattr(self, "_is_dual_index"):
            self._is_dual_index = {} 
            for lane in self.lanes:
                if self.index2cycles[lane] == "0": 
                    self._is_dual_index[lane] = False
                else:
                    self._is_dual_index[lane] = True 
        return self._is_dual_index

    @property
    def index1cycles(self):
        if not hasattr(self, "_index1cycles"):
            self._index1cycles = {} 
            for lane in self.lanes:
                self._index1cycles[lane] = self.get_index1cycles(lane)
        return self._index1cycles

    @property
    def index2cycles(self):
        if not hasattr(self, "_index2cycles"):
            self._index2cycles = {} 
            for lane in self.lanes:
                self._index2cycles[lane] = self.get_index2cycles(lane)
        return self._index2cycles

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
    def readset_file(self):
        if not hasattr(self, "_readset_file"):
            if self.args.readsets:
                self._readset_file = os.path.realpath(self.args.readsets.name)
            else:
                _raise(SanitycheckError("Error: missing '-r/--readsets' argument !"))
        return self._readset_file

    @property
    def lane_number(self):
        if self.args.lane_number:
            return self.args.lane_number
        else:
            return None

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
                lane_instrument = self.get_instrument(lane)
                if self._instrument and self._instrument != lane_instrument:
                    _raise(SanitycheckError("One run (\"" + self.run_id + "\") cannot be shared by different instruments !! (\"" + lane_instrument + "\" vs. \"" +  self._instrument + "\")"))
                else:
                    self._instrument = lane_instrument
        return self._instrument

    @property
    def flowcell_position(self):
        if not hasattr(self, "_flowcell_position"):
            self._flowcell_position = ""
            for lane in self.lanes:
                lane_flowcell_position = self.get_flowcell_position(lane)
                if self._flowcell_position and self._flowcell_position != lane_flowcell_position:
                    _raise(SanitycheckError("One run (\"" + self.run_id + "\") cannot be placed on multiple flowcell positions !! (\"" + lane_flowcell_position + "\" vs. \"" +  self._flowcell_position + "\")"))
                else:
                    self._flowcell_position = lane_flowcell_position
        return self._flowcell_position

    @property
    def run_number(self):
        if not hasattr(self, "_run_number"):
            run_number = ""
            for lane in self.lanes:
                lane_run_number = self.get_run_number(lane)
                if run_number and run_number != lane_run_number:
                    _raise(SanitycheckError("One run (\"" + self.run_id + "\") cannot be defined by more than one run counter !! (\"" + lane_run_number + "\" vs. \"" +  run_number + "\")"))
                else:
                    run_number = lane_run_number
            self._run_number = run_number
        return self._run_number

    @property
    def sequencer_run_id(self):
        if not hasattr(self, "_sequencer_run_id"):
            sequencer_run_id = "" 
            for lane in self.lanes:
                lane_sequencer_run_id = self.get_sequencer_run_id(lane)
                if sequencer_run_id and sequencer_run_id != lane_sequencer_run_id:
                    _raise(SanitycheckError("Sequencer Run ID conflct (\"" + lane_sequencer_run_id + "\" vs. \"" +  sequencer_run_id + "\")"))
                else:
                    sequencer_run_id = lane_sequencer_run_id
            self._sequencer_run_id = sequencer_run_id
        return self._sequencer_run_id

    @property
    def seqtype(self):
        if not hasattr(self, "_seqtype"):
            self._seqtype = self.get_seqtype()
        return self._seqtype

    @property
    def seq_category(self):
        if not hasattr(self, "_seq_category"):
            self._seq_category = self.get_seq_category()
        return self._seq_category

    @property
    def year(self):
        """
        Get year of the from sample sheet
        """
        if not hasattr(self, "_year"):
            dates = set([date for date in list(set([line['Start Date'] for line in csv.DictReader(open(self.readset_file, 'r'), delimiter='\t', quotechar='"')]))])
            if len(list(dates)) > 1:
                _raise(SanitycheckError("More than one date were found in the sample sheet for the run \"" + self.run_id + "\""))
            else:
                self._year = list(dates)[0].split("-")[0]
        return self._year

    @property
    def date(self):
        """
        Get whole date of the run from sample sheet
        """
        if not hasattr(self, "_date"):
            dates = set([date for date in list(set([line['Start Date'] for line in csv.DictReader(open(self.readset_file, 'r'), delimiter='\t', quotechar='"')]))])
            if len(list(dates)) > 1:
                _raise(SanitycheckError("More than one date were found in the sample sheet for the run \"" + self.run_id + "\""))
            else:
                date = list(dates)[0].split("-")
                self._date = date[0][-2:] + date[1] + date[2]
        return self._date

    @property
    def report_hash(self):
        if not hasattr(self, "_report_hash"):
            self._report_hash = {}
            for lane in self.lanes:
                self._report_hash[lane] = {
                    "version" : "1.0",
                    "run" : self.run_id,
                    "instrument" : self.instrument,
                    "flowcell" : self.flowcell_id,
                    "lane" : lane,
                    "seqtype" : self.seqtype,
                    "sequencing_method" : "PAIRED_END" if self.is_paired_end[lane] else "SINGLE_END",
                    "steps" : [],
                    "barcodes" : {}
                }
        return self._report_hash

    @property
    def report_inputs(self):
        if not hasattr(self, "_report_inputs"):
            self._report_inputs = {}
            for lane in self.lanes:
                self._report_inputs[lane] = {
                    'basecall': "",
                    'index': "",
                    'fastq': "",
                    'fastqc': {},
                    'qc_graphs' : {},
                    'blast' : {},
                    'picard_mark_duplicates' : {},
                    'metrics' : {}
                }
        return self._report_inputs

    @property
    def report_dir(self):
        if not hasattr(self, "_report_dir"):
            self._report_dir = {}
            for lane in self.lanes:
                self._report_dir[lane] = os.path.join(self.output_dir, "report")
        return self._report_dir

    @property
    def run_validation_report_json(self):
        if not hasattr(self, "_run_validation_report_json"):
            self._run_validation_report_json = {}
            for lane in self.lanes:
                self._run_validation_report_json[lane] = os.path.join(self.report_dir[lane], self.run_id + "." + lane + ".run_validation_report.json")
        return self._run_validation_report_json    

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
                region_name = readset.name + "_" + readset.sample_number + "_L00" + lane

                lane_jobs.append(
                    concat_jobs([
                        bash.mkdir(output_dir),
                        bvatools.readsqc(
                            readset.fastq1,
                            readset.fastq2,
                            "FASTQ",
                            region_name,
                            output_dir
                    )],
                    name="qc_graphs." + readset.name + ".qc." + self.run_id + "." + lane,
                    samples=[readset.sample]
                ))
                self.report_inputs[lane]['qc_graphs'][readset.name] = os.path.join(output_dir, "mpsQC_" + region_name + "_stats.xml")

                if not self.no_index_fastq:
                    if readset.index_fastq1:
                        # Also process the fastq of the indexes
                        lane_jobs.append(
                            concat_jobs(
                                [
                                    bash.mkdir(output_dir + "_index"),
                                    bvatools.readsqc(
                                        readset.index_fastq1,
                                        readset.index_fastq2 if self.is_dual_index[lane] else None,
                                        "FASTQ",
                                        region_name,
                                        output_dir + "_index"
                                    )
                                ],
                                name="qc_graphs." + readset.name + ".qc_index." + self.run_id + "." + lane,
                                samples=[readset.sample]
                            )
                        )

            self.add_to_report_hash("qc_graphs", lane, lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)

        if hasattr(self.args, 'type') and self.args.type == 't7':
            return jobs
        else:
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
                        os.path.join(os.path.dirname(readset.fastq1), "fastqc.R1", re.sub(".fastq.gz", "_fastqc", os.path.basename(readset.fastq1)), re.sub(".fastq.gz", "_fastqc", "fastqc_data.txt"))
                    ]
                }
                if readset.run_type == "PAIRED_END":
                    input_dict[readset.fastq2] = [
                        os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.fastq2))),
                        os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.fastq2))),
                        os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc", os.path.basename(readset.fastq2)), re.sub(".fastq.gz", "_fastqc", "fastqc_data.txt"))
                    ]
                if not self.no_index_fastq:
                    if readset.index_fastq1:
                        input_dict[readset.index_fastq1] = [
                            os.path.join(os.path.dirname(readset.index_fastq1), "fastqc.Barcodes", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.index_fastq1))),
                            os.path.join(os.path.dirname(readset.index_fastq1), "fastqc.Barcodes", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.index_fastq1)))
                        ]
                for input1, outputs in input_dict.items():
                    unzip = False
                    if input1 == readset.fastq1:
                        job_suffix = "R1."
                        input2 = None
                        unzip = True
                        if not len(self.readsets[lane]) == 1:
                            self.report_inputs[lane]['fastqc'][readset.name] = outputs[2]
                    elif input1 == readset.fastq2:
                        job_suffix = "R2."
                        input2 = None
                        unzip=True
                    elif input1 == readset.index_fastq1:
                        job_suffix = "Barcodes."
                        input2 = readset.index_fastq2 if self.is_dual_index[lane] else None
                    lane_jobs.append(
                        concat_jobs([
                            bash.mkdir(
                                os.path.dirname(outputs[0]),
                                remove=True
                            ),
                            fastqc.fastqc(
                                input1,
                                input2,
                                outputs,
                                adapter_file=None,
                                extract=unzip,
                                use_tmp=True
                        )],
                        name="fastqc." + readset.name + "_" + job_suffix + "." + self.run_id + "." + lane,
                        samples=[readset.sample]
                    ))

            self.add_to_report_hash("fastqc", lane, lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)

        if hasattr(self.args, 'type') and self.args.type == 't7':
            return jobs
        else:
            return self.throttle_jobs(jobs)

    def fastp(self):
        """
        Generate basic QC metrics.
        """
        jobs = []

        for lane in self.lanes:
            lane_jobs = []
            for readset in self.readsets[lane]:
                output_json_path = os.path.join(os.path.dirname(readset.fastq1), "fastp", readset.name + ".fastp.json")
                output_html_path = os.path.join(os.path.dirname(readset.fastq1), "fastp", readset.name + ".fastp.html")
                lane_jobs.append(
                    concat_jobs([
                        bash.mkdir(os.path.dirname(output_json_path), remove=True),
                        fastp.fastp_basic_qc(readset.fastq1, readset.fastq2, output_json_path, output_html_path),
                    ],
                    name="fastp." + readset.name + "." + self.run_id + "." + lane,
                    samples=[readset.sample]
                    )
                )
            self.add_to_report_hash("fastp", lane, lane_jobs)
            jobs.extend(lane_jobs)
        return jobs

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

        for lane in self.lanes:
            lane_jobs = []

            if is_nb_blast_per_lane:
                nb_blast_to_do = int(nb_blast_to_do) // len(self.readsets[lane])

            nb_blast_to_do = max(1, nb_blast_to_do)

            for readset in self.readsets[lane]:
                output_prefix = os.path.join(
                    self.output_dir,
                    "Unaligned." + lane,
                    "Blast_sample",
                    readset.name + "_" + readset.sample_number + "_L00" + lane
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
                    # Because runBlast.sh ends up creating a symlink, here
                    #  we make sure to remove the link before it's created.
                    # We also add the R2.fastq as a parameter to runBlast.sh
                    command = "rm -f " + output + " && " + command + readset.fastq2
                    fasta_file = output_prefix + ".R1R2.subSampled_{nb_blast_to_do}.fasta".format(
                        nb_blast_to_do=nb_blast_to_do
                    )
                    result_file = output_prefix + ".R1R2.subSampled_{nb_blast_to_do}.blastres".format(
                        nb_blast_to_do=nb_blast_to_do
                    )
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
                    if readset.fastq2:
                        rrna_output = output_prefix + ".R1R2.subSampled_{nb_blast_to_do}.rrna".format(
                            nb_blast_to_do=nb_blast_to_do
                        )
                    else:
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

            self.add_to_report_hash("blast", lane, lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)

        if hasattr(self.args, 'type') and self.args.type == 't7':
            return jobs
        else:
            return self.throttle_jobs(jobs)

    def align(self):
        """
        Align the reads from the fastq file, sort the resulting .bam and create an index
        of that .bam.

        An basic aligment is performed on a sample when the `SampleRef` field of the
        MGI sample sheet match one of the regexp in the configuration file and the
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

            self.add_to_report_hash("align", lane, lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)

        if hasattr(self.args, 'type') and self.args.type == 't7':
            return jobs
        else:
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
                self.report_inputs[lane]['picard_mark_duplicates'][readset.name] = metrics_file

                job = picard.mark_duplicates(
                    [input],
                    output,
                    metrics_file
                )
                job.name = "picard_mark_duplicates." + readset.name + ".dup." + self.run_id + "." + lane
                job.samples = [readset.sample]
                lane_jobs.append(job)

            self.add_to_report_hash("picard_mark_duplicates", lane, lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)

            jobs.extend(lane_jobs)

        if hasattr(self.args, 'type') and self.args.type == 't7':
            return jobs
        else:
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
                    self.report_inputs[lane]['metrics'][readset.name] = [
                        readset.bam + '.metrics.alignment_summary_metrics',
                        readset.bam + '.metrics.insert_size_metrics',
                        os.path.join(os.path.dirname(readset.bam), readset.sample.name + "." + readset.library + '.rnaseqc.sorted.dup.metrics.tsv'),
                        readset.bam + ".metrics.verifyBamId.tsv",
                        readset.bam + '.metrics.rRNA.tsv'
                    ]
                else:
                    self.report_inputs[lane]['metrics'][readset.name] = [
                        readset.bam + '.metrics.alignment_summary_metrics',
                        readset.bam + '.metrics.insert_size_metrics',
                        readset.bam + ".metrics.verifyBamId.tsv",
                        readset.bam + ".metrics.targetCoverage.txt"
                    ]

            self.add_to_report_hash("metrics", lane, lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)

        if hasattr(self.args, 'type') and self.args.type == 't7':
            return jobs
        else:
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
                    bash.md5sum(
                        readset.fastq1,
                        readset.fastq1 + ".md5",
                        binary=True
                    )
                ]

                # Second read in paired-end run
                if readset.fastq2:
                    current_jobs.append(
                        bash.md5sum(
                            readset.fastq2,
                            readset.fastq2 + ".md5",
                            binary=True
                        )
                    )

                # Alignment files
                if readset.bam:
                    current_jobs.append(
                        bash.md5sum(
                            readset.bam + ".bam",
                            readset.bam + ".bam.md5",
                            binary=True
                        )
                    )
                    current_jobs.append(
                        bash.md5sum(
                            readset.bam + ".bai",
                            readset.bam + ".bai.md5",
                            binary=True
                        )
                    )

                job = concat_jobs(
                    current_jobs,
                    name="md5." + readset.name + ".md5." + self.run_id + "." + lane,
                    samples=[readset.sample]
                )

                lane_jobs.append(job)

            if config.param('md5', 'one_job', required=False, type="boolean"):
                lane_job = concat_jobs(
                    lane_jobs,
                    name="md5." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
                self.add_copy_job_inputs([lane_job], lane)
                jobs.append(lane_job)
            else:
                self.add_copy_job_inputs(lane_jobs, lane)
                jobs.extend(lane_jobs)

        return jobs

    def report(self):
        """
        Generate a JSON file reporting the whole pipeline
        """
        jobs = []

        for lane in self.lanes:
            lane_jobs = []

            report_dir = os.path.join(self.output_dir, "report")
            sample_report_dir  = os.path.join(report_dir, "sample_json")

            run_validation_inputs = []

            # Add barcodes info to the report_hash
            self.report_hash[lane]["barcodes"] = dict([(readset.name, readset.indexes) for readset in self.readsets[lane]])

            self.generate_lane_json_report_file(lane)

            general_information_file = os.path.join(self.output_dir, self.run_id + "." + lane + ".general_information.json")
            if not os.path.exists(os.path.dirname(general_information_file)):
                os.makedirs(os.path.dirname(general_information_file))
            with open(general_information_file, 'w') as out_json:
                json.dump(self.report_hash[lane], out_json, indent=4)

            # metrics to JSON
            for step in self.step_list:
                report_step_jobs = []
                if step.name in self.report_inputs[lane].keys() and self.report_inputs[lane][step.name]:
                   pipeline = self.args.type if hasattr(self.args, 'type') else 'illumina'
                   if type(self.report_inputs[lane][step.name]) is str:
                       report_step_jobs.append(
                           concat_jobs(
                               [
                                   tools.run_processing_metrics_to_json(
                                       self.run_validation_report_json[lane],
                                       step.name,
                                       pipeline,
                                       self.report_inputs[lane][step.name]
                                   ),
                                   bash.touch(os.path.join(self.job_output_dir, "checkpoint", step.name + "." + self.run_id + "." + lane + ".reportUpdated"))
                               ],
                               name="report." + step.name + "." + self.run_id + "." + lane,
                               samples=self.samples[lane]
                           )
                       )
                   elif type(self.report_inputs[lane][step.name]) is dict:
                       for readset in self.readsets[lane]:
                           report_step_jobs.append(
                               concat_jobs(
                                   [
                                       tools.run_processing_metrics_to_json(
                                           self.run_validation_report_json[lane],
                                           step.name,
                                           pipeline,
                                           self.report_inputs[lane][step.name][readset.name],
                                           readset.name
                                       ),
                                       bash.touch(os.path.join(self.job_output_dir, "checkpoint", step.name + "." + readset.name + "." + self.run_id + "." + lane + ".reportUpdated"))
                                   ],
                                   name="report." + step.name + "." + readset.name + "." + self.run_id + "." + lane,
                                   samples=self.samples[lane]
                               )
                           )
                   else:
                       _raise(SanitycheckError("Unknown metrics type : " + type(self.report_inputs[lane][step.name]) + " for step '" + step.name + "'"))

                   lane_jobs.extend(report_step_jobs)

                # checkpoint file
                step_checkpoint_file = os.path.join(self.job_output_dir, "checkpoint", step.name + "." + self.run_id + "." + lane + ".done")
                if report_step_jobs:
                    step_checkpoint_job_dependencies = [output for job in report_step_jobs for output in job.output_files]
                else:
                    step_checkpoint_job_dependencies = [output for job in step.jobs for output in job.output_files]
                if step_checkpoint_job_dependencies:
                    lane_jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(os.path.dirname(step_checkpoint_file)),
                                bash.touch(step_checkpoint_file)
                            ],
                            name="checkpoint." + step.name + "." + self.run_id + "." + lane,
                            input_dependency=step_checkpoint_job_dependencies,
                            output_dependency=[step_checkpoint_file],
                            samples=self.samples[lane]
                        )
                    )

            self.add_copy_job_inputs(lane_jobs, lane)

            jobs.extend(lane_jobs)

        return jobs

    def copy(self):
        """
        Copy the whole processing foler to where they can be serve or loaded into a LIMS
        """
        jobs_to_concat = []

        if hasattr(self.args, 'type'):
            full_destination_folder = os.path.join(
                config.param("copy", "destination_folder", type="dirpath"),
                self.seq_category,
                self.year,
                self.date + "_" + self.instrument + "_" + self.run_number + "_" + self.flowcell_position + self.flowcell_id + "_" + self.sequencer_run_id + "-" + self.seqtype
            )
        else:
            full_destination_folder = os.path.join(
                config.param("copy", "destination_folder", type="dirpath"),
                self.seq_category,
                self.year,
                os.path.basename(self.run_dir) + "-" + self.seqtype
            )

        jobs_to_concat.append(
            bash.mkdir(full_destination_folder)
        )

        jobs_to_concat.append(
            bash.cp(self.readset_file, full_destination_folder)
        )

        for lane in self.lanes:

            inputs = self.copy_job_inputs[lane]
            copy_output = os.path.join(
                full_destination_folder,
                "copyCompleted." + lane + ".out"
            )
            report_transfer_output = os.path.join(
                full_destination_folder,
                "reportTransferCompleted." + lane + ".out"
            )

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

            copy_command_output_folder = config.param('copy', 'copy_command', required=False).format(
                exclusion_clauses="\\\n".join(
                    [" --exclude '" + excludedfile.replace(self.output_dir + os.sep, "") + "'" for excludedfile in excluded_files]),
                lane_number=lane,
                run_id=self.run_id,
                source=os.path.join(self.output_dir),
                destination=full_destination_folder
            )

            jobs_to_concat.append(concat_jobs(
                [
                    Job(
                        command=copy_command_output_folder
                    ),
                    bash.touch(copy_output)
                ],
                input_dependency=inputs,
                output_dependency=[copy_output]
            ))

        job = concat_jobs(
            jobs_to_concat,
            name="copy." + self.run_id,
            samples=self.samples[lane]
        )
        return [job]

    def final_notification(self):
        """
        Writes a simple '.done' file when all pipeline is done processing
        """
        jobs = []

        inputs = []
        for lane in self.lanes:
            if hasattr(self.args, 'type'):
                full_destination_folder = os.path.join(
                    config.param("copy", "destination_folder", type="dirpath"),
                    self.seq_category,
                    self.year,
                    self.date + "_" + self.instrument + "_" + self.run_number + "_" + self.flowcell_position + self.flowcell_id + "_" + self.sequencer_run_id + "-" + self.seqtype
                )
            else:
                full_destination_folder = os.path.join(
                    config.param("copy", "destination_folder", type="dirpath"),
                    self.seq_category,
                    self.year,
                    os.path.basename(self.run_dir) + "-" + self.seqtype
                )

            inputs.append(os.path.join(full_destination_folder, "copyCompleted." + lane + ".out"))

        notification_job = bash.touch(os.path.join(self.output_dir,  self.run_id + "_" + self.flowcell_id + "_processing.done"))
        notification_job.input_files = inputs
        notification_job.output_files = [os.path.join(self.output_dir,  self.run_id + "_" + self.flowcell_id + "_processing.done")]
        notification_job.name = "final_notification." + self.run_id
        notification_job.samples = self.samples[lane]

        return [notification_job]

    #
    # Utility methods
    #

    def add_copy_job_inputs(self, jobs, lane):
        for job in jobs:
            # we first remove dependencies of the current job, since we will have a dependency on that job
            self.copy_job_inputs[lane] = [item for item in self.copy_job_inputs[lane] if item not in job.input_files]
            self.copy_job_inputs[lane].extend(job.output_files)

    def add_to_report_hash(self, step_name, lane, jobs=[]):
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
        self.report_hash[lane]["steps"].append(report_hash)

    def generate_lane_json_report_file(self, lane):
        """
        Builds a JSON object containing :
           - general information about the run
           - list and description of the pipeline jobs
           - pipeline metrics organized by step
        Dump the JSON object into a report file
        """

        self.report_hash[lane]["total_pf_clusters"] = None
        self.report_hash[lane]["spread"] = None
        self.report_hash[lane]["run_validation"] = []
        for readset in self.readsets[lane]:
            self.report_hash[lane]["run_validation"].append(
                {
                    "project": readset.project_id,
                    "sample": readset.name,
                    "index": {
                        "Barcode": readset.index_name,
                        "Barcode sequence": ','.join([readset_index['BARCODE_SEQUENCE'] for readset_index in readset.indexes]),
                        "% on index in lane": None,
                        "% of the lane": None,
                        "% Perfect barcode": None,
                        "% One mismatch barcode": None,
                        "PF Clusters": None,
                        "Yield (bases)": None,
                        "Mean Quality Score": None,
                        "% >= Q30 bases": None
                    },
                    "qc": {
                        "avgQual": None,
                        "duplicateRate": None
                    },
                    "sample_tag": None,
                    "blast": {
                        "1st_hit": None,
                        "2nd_hit": None,
                        "3rd_hit": None
                    },
                    "alignment": {
                        "chimeras": None,
                        "average_aligned_insert_size": None,
                        "reported_sex": readset.gender,
                        "pf_read_alignment_rate": None,
                        "Freemix": None,
                        "inferred_sex": None,
                        "adapter_dimers": None,
                        "mean_coverage": None,
                        "aligned_dup_rate": None,
                        "sex_concordance": None
                    }
                }
            )

        report_dir = os.path.join(self.output_dir, "report")
        if not os.path.exists(os.path.dirname(self.run_validation_report_json[lane])):
             os.makedirs(os.path.dirname(self.run_validation_report_json[lane]))
        with open(self.run_validation_report_json[lane], 'w') as out_json:
            json.dump(self.report_hash[lane], out_json, indent=4)

    def get_seqtype(self):
        """
        Determine which kind of sequencing (iseq, miseq, novaseq, hiseqx, hiseq4000, hiseq2500 or dnbseqg400) was performed,
        depending on the instrument used for the run
        """

        instrument = self.instrument
        instrument_file = config.param('DEFAULT', 'instrument_list_file', type='filepath', required=False)
        if not (instrument_file and os.path.isfile(instrument_file)):
            instrument_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'instrument_list.csv')

        return subprocess.check_output("grep -m1 '"+instrument+"' %s | awk -F',' '{print $3}'" % instrument_file, shell=True, text=True).strip()

    def get_seq_category(self):
        if "hiseq" in self.seqtype:
            return "hiseq"
        else:
            return self.seqtype

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
                number_task_by_job = int(math.ceil(len(current_jobs) / float(max_jobs_per_step)))
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
 
