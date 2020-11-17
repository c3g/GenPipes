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
import csv
import collections
import json

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs

from bfx.readset import parse_mgi_raw_readset_files
from bfx import bvatools
from bfx import picard
from bfx import fastqc
from bfx import tools
from bfx import run_processing_tools
from bfx import bash_cmd as bash

from pipelines import common

log = logging.getLogger(__name__)

class MGIRunProcessing(common.MUGQICPipeline):
    """
    MGI Run Processing Pipeline
    ================================

    The standard MUGQIC MGI Run Processing pipeline uses fastq files produced
    by the sequencer, then does demultiplexing. Finally, the
    pipeline runs some QCs on the raw data, on the fastq and on the alignment.

    Sample Sheets
    -------------

    The pipeline uses one input sample sheet.
    CURRENTLY BASED ON MGI RUN PROCESSING GOOGLE SHEET (https://docs.google.com/spreadsheets/d/1Jk11bQUJdqVg37gfn7ndfk-g9ke96tsjCfAsf3r1xdA)
    A csv file having the following columns :

    - Sample
    - Readset
    - Library
    - Project
    - Project ID
    - Protocol
    - Index
    - Pool ID
    - Run ID
    - Flowcell ID
    - Lane
    - Sequencer
    - Sequencer ID

    Example:
        Sample,Readset,Library,Project,Project ID,Protocol,Index,PoolID,RunID,FlowcellID,Lane,Sequencer,SequencerID
        LSPQ_Viral_Culture_dil_10-1_10cycles,LSPQ_Viral_Culture_dil_10-1_10cycles_PROD_000034-A01,PROD_000034-A01,LSPQ,,CleanPlex_MGI,1,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        LSPQ_Viral_Culture_dil_10-2_10cycles,LSPQ_Viral_Culture_dil_10-2_10cycles_PROD_000034-B01,PROD_000034-B01,LSPQ,,CleanPlex_MGI,2,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        LSPQ_Viral_Culture_dil_10-3_10cycles,LSPQ_Viral_Culture_dil_10-3_10cycles_PROD_000034-C01,PROD_000034-C01,LSPQ,,CleanPlex_MGI,3,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        LSPQ_Nasal_Swab_Neg_ctl_10cycles,LSPQ_Nasal_Swab_Neg_ctl_10cycles_PROD_000034-A02,PROD_000034-A02,LSPQ,,CleanPlex_MGI,25,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        LSPQ_Viral_Culture_dil_10-4_13cycles,LSPQ_Viral_Culture_dil_10-4_13cycles_PROD_000034-A03,PROD_000034-A03,LSPQ,,CleanPlex_MGI,28,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        LSPQ_Viral_Culture_dil_10-5_13cycles,LSPQ_Viral_Culture_dil_10-5_13cycles_PROD_000034-B03,PROD_000034-B03,LSPQ,,CleanPlex_MGI,29,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        L00241026_dil_10-1_13cycles,L00241026_dil_10-1_13cycles_PROD_000034-E03,PROD_000034-E03,LSPQ,,CleanPlex_MGI,33,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        L00241026_dil_10-2_13cycles,L00241026_dil_10-2_13cycles_PROD_000034-F03,PROD_000034-F03,LSPQ,,CleanPlex_MGI,34,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        Arctic_RT_reaction_13cycles,Arctic_RT_reaction_13cycles_PROD_000034-B04,PROD_000034-B04,LSPQ,,CleanPlex_MGI,4,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01

    """

    def __init__(self, protocol=None):
        self._protocol=protocol
        self.copy_job_inputs = []
        self.argparser.add_argument("-d", "--run", help="Run directory (mandatory)", required=True, dest="run_dir")
        self.argparser.add_argument("--lane", help="Lane number (mandatory)", type=int, required=True, dest="lane_number")
        self.argparser.add_argument("-r", "--readsets", help="Sample sheet for the MGI run to process (mandatory)", type=file, required=True)
        self.argparser.add_argument("--raw-fastq", help="Raw fastq from the sequencer, with NO DEMULTIPLEXING performed, will be expected by the pipeline", action="store_true", required=False, dest="raw_fastq")
        self.argparser.add_argument("-x", help="First index base to use for demultiplexing (inclusive). The index from the sample sheet will be adjusted according to that value.", type=int, required=False, dest="first_index")
        self.argparser.add_argument("-y", help="Last index base to use for demultiplexing (inclusive)", type=int, required=False, dest="last_index")
        self.argparser.add_argument("-m", help="Number of index mistmaches allowed for demultiplexing (default 1). Barcode collisions are always checked.", type=int, required=False, dest="number_of_mismatches")

        super(MGIRunProcessing, self).__init__(protocol)

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = self.load_readsets()
            if not self.is_demultiplexed and len(self.readsets) > 1:
                self.generate_mgi_lane_barcode_sheet()
                self.generate_mgi_lane_sample_sheet()
        return self._readsets

    @property
    def is_paired_end(self):
        if not hasattr(self, "_is_paired_end"):
            if self.get_read2cycles():
                self._is_paired_end = True
            else:
                self._is_paired_end = False
        return self._is_paired_end

    @property
    def is_dual_index(self):
        if not hasattr(self, "_is_dual_index"):
            if self.get_index2cycles() == "0":
                self._is_dual_index = False
            else:
                self._is_dual_index = True
        return self._is_dual_index

    @property
    def is_demultiplexed(self):
        if not hasattr(self, "_is_demultiplexed"):
            if self.args.raw_fastq:
                self._is_demultiplexed = False
            else:
                self._is_demultiplexed = True
        return self._is_demultiplexed

    @property
    def mask(self):
        if not hasattr(self, "_mask"):
            self._mask = self.get_mask()
        return self._mask 

    @property
    def merge_undetermined(self):
        if not hasattr(self, "_merge_undetermined"):
            self._merge_undetermined = False
            # If only one library on the lane
            if len(self.readsets) == 1:
                self._merge_undetermined = True
            if config.param('fastq', 'merge_undetermined', required=False, type='boolean'):
                self._merge_undetermined = config.param('fastq', 'merge_undetermined')
        return self._merge_undetermined

    @property
    def instrument(self):
        if not hasattr(self, "_instrument"):
            self._instrument = self.get_instrument()
        return self._instrument

    @property
    def seqtype(self):
        if not hasattr(self, "_seqtype"):
            self._seqtype = self.get_seqtype()
        return self._seqtype

    @property
    def year(self):
        """
        """
        if not hasattr(self, "_year"):
            self._year = "2020"
        return self._year 

    @property
    def run_id(self):
        """
        The run id from the readset objects.
        """
        if not hasattr(self, "_run_id"):
            runs = set([readset.run for readset in self.readsets])
            if len(list(runs)) > 1:
                _raise(SanitycheckError("Error: more than one run were parsed in the sample sheet... " + runs))
            else:
                self._run_id = list(runs)[0]
        return self._run_id

    @property
    def run_dir(self):
        if self.args.run_dir:
            return self.args.run_dir
        else:
            _raise(SanitycheckError("Error: missing '-d/--run_dir' option!"))

    @property
    def bioinfo_file(self):
        if not hasattr(self, "_bioinfo_file"):
            if os.path.exists(os.path.join(self.output_dir, "Unaligned." + str(self.lane_number), "raw_fastq", "BioInfo.csv")):
                self._bioinfo_file = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number), "raw_fastq", "BioInfo.csv")
            else:
                rundir = self.run_dir
                if "(no_barcode)" in rundir:
                    rundir = rundir.split("(")[0]
                self._bioinfo_file = os.path.join(rundir, "L0" + str(self.lane_number), "BioInfo.csv")
        return self._bioinfo_file

    @property
    def flowcell_id(self):
        """
        The flow cell ID from the run folder
        """
        if not hasattr(self, "_flowcell_id"):
            self._flowcell_id = os.path.basename(self.run_dir.rstrip('/'))
            if "(" in self._flowcell_id:
                 self._flowcell_id = self._flowcell_id.split("(")[0]
        return self._flowcell_id

    @property
    def lane_number(self):
        if self.args.lane_number:
            return self.args.lane_number
        else:
            _raise(SanitycheckError("Error: missing '--lane' option!"))

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
    def index_per_readset(self):
        return self._index_per_readset

    @property
    def readset_file(self):
        if self.args.readsets:
            return self.args.readsets.name
        else:
            _raise(SanitycheckError("Error: missing '-r/--readsets' argument !"))

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
                "barcodes" : {}
            }
        return self._report_hash
    
    @property
    def report_inputs(self):
        if not hasattr(self, "_report_inputs"):
            self._report_inputs = {
                'qc' : {},
                'blast' : {},
                'mark_dup' : {},
                'align' : {}
            }
        return self._report_inputs

    def index(self):
        """
        First copy all the files of the lane from the sequencer deposit folder to the processing folder,
        into "raw_fastq".
        Then, if MGI adapters are used (i.e. demultiplexing already done on sequencer)) then
            rename the fastq files and move them from "raw_fastq" to "Unaligned.LANE" folder,
            and pipeline will SKIP demultiplexing (i.e. next step, fastq).
        Else, for all non-MGI adapters (i.e. demultiplexing still remains to be done) then prepare
            prepare data for demultiplexing (i.e. next step, fastq) by extracting adapters from R2 fastq,
            creating I1 and I2 fastq files for adapter i7 and i5 respectively (or only I1 form i7 in case
            of single index sequencing). This will also create the barcode tsv file used to launch
            fastq-multx to demultiplex the reads in fastq step.
        *TO DO* - in both cases, retrieve index stats from the sequencer output files to build an index report
        """

        jobs = []

        # First we copy all the lane folder into the raw_fastq folder
        unaligned_dir = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number))
        raw_fastq_dir = os.path.join(unaligned_dir, 'raw_fastq')
        copy_done_file = os.path.join(raw_fastq_dir, "copy_done.Success")
        copy_job_output_dependency = [
            os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + ".summaryReport.html"),
            os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + ".heatmapReport.html"),
            os.path.join(raw_fastq_dir, "summaryTable.csv")
        ]
        rename_job = None

        # If demultiplexing were done on the sequencer i.e. MGI adpaters used...
        if self.is_demultiplexed:
            fastq_jobs = []

            # ...then do the necessary moves and renamings, from the raw_fastq folder to the Unaligned folders
            for readset in self.readsets:
                output_dir = os.path.join(self.output_dir, "Unaligned." + readset.lane, 'Project_' + readset.project, 'Sample_' + readset.name)
                m = re.search("\D+(?P<index>\d+)", readset.index_name)
                if m:
                    index_number = m.group('index').lstrip("0")
                else:
                    _raise(SanitycheckError("MGI Index error : unable to parse barcode # in : " + readset.index_name))

                copy_job_output_dependency.extend([
                    os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + index_number + "_1.fq.gz"),
                    os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + index_number + "_1.fq.fqStat.txt"),
                    os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + index_number + ".report.html")
                ])
                fastq_job_input_dependency.extend([
                    os.path.join(self.run_dir, "L0" + str(self.lane_number), readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_1.fq.gz"),
                    os.path.join(self.run_dir, "L0" + str(self.lane_number), readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_1.fq.fqStat.txt"),
                    os.path.join(self.run_dir, "L0" + str(self.lane_number), readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + ".report.html")
                ])
                fastq_jobs.append(
                    concat_jobs([
                        bash.mkdir(output_dir),
                        bash.cp(
                            os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + index_number + "_1.fq.gz"),
                            readset.fastq1
                        ),
                        bash.cp(
                            os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + index_number + "_1.fq.fqStat.txt"),
                            re.sub("gz", "fqStat.txt", readset.fastq1)
                        ),
                        bash.cp(
                            os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + index_number + ".report.html"),
                            re.sub("_R1_001.fastq.gz", ".report.html", readset.fastq1)
                        ),
                        bash.ln(
                            re.sub("_R1_001.fastq.gz", ".report.html", readset.fastq1),
                            os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + ".report.html")
                        )
                    ])
                )

                if readset.run_type == "PAIRED_END":
                    copy_job_output_dependency.extend([
                        os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + index_number + "_2.fq.gz"),
                        os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + index_number + "_2.fq.fqStat.txt")
                    ])
                    fastq_job_input_dependency.extend([
                        os.path.join(self.run_dir, "L0" + str(self.lane_number), readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_2.fq.gz"),
                        os.path.join(self.run_dir, "L0" + str(self.lane_number), readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_2.fq.fqStat.txt")
                    ])
                    fastq_jobs.append(
                        concat_jobs([
                            bash.cp(
                                os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + index_number + "_2.fq.gz"),
                                readset.fastq2,
                            ),
                            bash.cp(
                                os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + index_number + "_2.fq.fqStat.txt"),
                                re.sub("gz", "fqStat.txt", readset.fastq2),
                            ),
                            bash.ln(
                                re.sub("gz", "fqStat.txt", readset.fastq2),
                                os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_2.fq.fqStat.txt")
                            )
                        ])
                    )
            fastq_job = concat_jobs(fastq_jobs)
            fastq_job.input_files = fastq_job_input_dependency + [copy_done_file]
            fastq_job.name = "index.copy_rename_raw." + self.run_id + "." + str(self.lane_number)
            fastq_job.samples = self.samples

        elif "(no_barcode)" in self.run_dir:
            if len(self.readsets) > 1:
                err_msg = "LANE SETTING ERROR :\n"
                err_msg += "Unable to demultiplex " + str(len(self.readsets)) + " samples : No barcode in fastq files...\n(in "
                err_msg += self.run_dir + ")"
                _raise(SanitycheckError(err_msg))
            readset = self.readsets[0]
            rename_job = concat_jobs(
                [
                    bash.mkdir(os.path.dirname(readset.fastq1)),
                    bash.cp(
                        os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_1.fq.gz"),
                        readset.fastq1,
                    ),
                    bash.cp(
                        os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_1.fq.fqStat.txt"),
                        re.sub("gz", "fqStat.txt", readset.fastq1),
                    ),
                    bash.cp(
                        os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read.report.html"),
                        re.sub("_R1_001.fastq.gz", "_read.report.html", readset.fastq1),
                    )
                ]
            )
            if readset.run_type == "PAIRED_END":
                rename_job = concat_jobs(
                    [
                        rename_job,
                        bash.cp(
                            os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz"),
                            readset.fastq2,
                        ),
                        bash.cp(
                            os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.fqStat.txt"),
                            re.sub("gz", "fqStat.txt", readset.fastq2),
                        )
                    ]
                )
            rename_job.output_files.append(copy_done_file)
            rename_job.name = "index.rename." + self.run_id + "." + str(self.lane_number)
            rename_job.samples = self.samples

        # If no demultiplexing were done on the sequencer i.e. no MGI adpater used... 
        else:
            # ...then extract read2, index1 and index2 from R2 fastq, and do some renamings before demultiplexing ('fastq' step)
            if len(self.readsets) == 1:
                readset = self.readsets[0]
                R1_fastq = readset.fastq1
                R2_fastq = readset.fastq2 if self.is_paired_end else None
                I1_fastq = readset.index_fastq1 
                I2_fastq = readset.index_fastq2 if self.is_dual_index else None
            else:
                R1_fastq = os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_R1.fastq.gz")
                R2_fastq = os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_R2.fastq.gz") if self.is_paired_end else None
                I1_fastq = os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_I1.fastq.gz")
                I2_fastq = os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_I2.fastq.gz") if self.is_dual_index else None

            # R1 fastq
            copy_job_output_dependency.append(R1_fastq)
            fastq_job_input_dependency.append(os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_1.fq.gz"))
            fastq_job = [concat_jobs(
                [
                    bash.mkdir(os.path.dirname(R1_fastq)),
                    bash.mv(
                        os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_1.fq.gz"),
                        R1_fastq
                    ),
                    bash.ln(
                        R1_fastq,
                        os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_1.fq.gz")
                    )
                ],
                input_dependency=[os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_1.fq.gz")],
                output_dependency=[R1_fastq],
                name="index.rename_R1_fastq." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            )]
            if self.is_paired_end:
                # R2 fastq
                copy_job_output_dependency.append(R2_fastq)
                fastq_job_input_dependency.append(os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz"))
                fastq_job.append(
                    pipe_jobs(
                        [
                            bash.cat(
                                os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz"),
                                None,
                                zip=True
                            ),
                            bash.cut(
                                None,
                                None,
                                options="-c 1-" + self.get_read2cycles()
                            ),
                            bash.gzip(
                                None,
                                R2_fastq
                            )
                        ],
                        input_dependency=[os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz")],
                        name="index.extract_R2_fastq." + self.run_id + "." + str(self.lane_number),
                        samples=self.samples
                    )
                )
            if self.is_dual_index:
                # I1 fastq
                copy_job_output_dependency.append(I1_fastq)
                fastq_job_input_dependency.append(os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz"))
                fastq_job.append(
                    pipe_jobs(
                        [
                            bash.cat(
                                os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz"),
                                None,
                                zip=True
                            ),
                            bash.paste(
                                None,
                                None,
                                options="-d '\\t' - - - -"
                            ),
                            bash.awk(
                                None,
                                None,
                                "-F'\\t' '{print $1 \"\\n\" substr($2,"+str(int(self.get_read2cycles())+int(self.get_index2cycles())+1)+","+self.get_index1cycles()+") \"\\n\" $3 \"\\n\" substr($4,"+str(int(self.get_read2cycles())+int(self.get_index2cycles())+1)+","+self.get_index1cycles()+") }'"
                            ),
                            bash.gzip(
                                None,
                                I1_fastq
                            )
                        ],
                        input_dependency=[os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz")],
                        name="index.extract_I1_fastq." + self.run_id + "." + str(self.lane_number),
                        samples=self.samples
                    )
                )
                # I2 fastq
                copy_job_output_dependency.append(I2_fastq)
                fastq_job_input_dependency.append(os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz"))
                fastq_job.append(
                    pipe_jobs(
                        [
                            bash.cat(
                                os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz"),
                                None,
                                zip=True
                            ),
                            bash.paste(
                                None,
                                None,
                                options="-d '\\t' - - - -"
                            ),
                            bash.awk(
                                None,
                                None,
                                "-F'\\t' '{print $1 \"\\n\" substr($2,"+str(int(self.get_read2cycles())+1)+","+self.get_index2cycles()+") \"\\n\" $3 \"\\n\" substr($4,"+str(int(self.get_read2cycles())+1)+","+self.get_index2cycles()+") }'"
                            ),
                            bash.gzip(
                                None,
                                I2_fastq
                            )
                        ],
                        input_dependency=[os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz")],
                        name="index.extract_I2_fastq." + self.run_id + "." + str(self.lane_number),
                        samples=self.samples
                    )
                )
            else:
                # I1 fastq
                copy_job_output_dependency.append(I1_fastq)
                fastq_job_input_dependency.append(os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz"))
                fastq_job.append(
                    pipe_jobs(
                        [
                            bash.cat(
                                os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz"),
                                None,
                                zip=True
                            ),
                            bash.paste(
                                None,
                                None,
                                options="-d '\\t' - - - -"
                            ),
                            bash.awk(
                                None,
                                None,
                                "-F'\\t' '{print $1 \"\\n\" substr($2,"+str(int(self.get_read2cycles())+1)+","+self.get_index1cycles()+") \"\\n\" $3 \"\\n\" substr($4,"+str(int(self.get_read2cycles())+1)+","+self.get_index1cycles()+") }'"
                            ),
                            bash.gzip(
                                None,
                                I1_fastq
                            )
                        ],
                        input_dependency=[os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz")],
                        name="index.extract_I1_fastq." + self.run_id + "." + str(self.lane_number),
                        samples=self.samples
                    )
                )

        # Job submission :
        # First we copy all the lane folder into the raw_fastq folder (unless already done...)
        if not os.path.exists(copy_done_file):
            copy_job = concat_jobs(
                [
                    bash.mkdir(raw_fastq_dir),
                    bash.cp(
                        os.path.join(self.run_dir, "L0" + str(self.lane_number), "."),
                        raw_fastq_dir,
                        recursive=True
                    ),
                    pipe_jobs([
                        bash.md5sum(
                            [re.sub(raw_fastq_dir, os.path.dirname(file_path), file_path) for file_path in copy_job_output_dependency],
                            None
                        ),
                        bash.awk(
                            None,
                            None,
                            "'sub( /\/.*\//,\"\",$2 )'"
                        ),
                        bash.awk(
                            None,
                            copy_done_file + ".md5",
                            "'{print $1,\""+raw_fastq_dir+"/\"$2}'"
                        )
                    ])
                ]
            )
            # Handle corner cases for BioInfo.csv
            if "(no_barcode)" in self.run_dir:
                copy_job = concat_jobs(
                    [
                        copy_job,
                        bash.cp(
                            self.bioinfo_file,
                            raw_fastq_dir
                        )
                    ]
                )
            copy_job = concat_jobs(
                [
                    copy_job,
                    pipe_jobs([
                        bash.md5sum(
                            self.bioinfo_file,
                            None
                        ),
                        bash.awk(
                             None,
                             None,
                             "'sub( /\/.*\//,\"\",$2 )'"
                        ),
                        bash.awk(
                            None,
                            copy_done_file + ".md5",
                            "'{print $1,\""+raw_fastq_dir+"/\"$2}'",
                            append=True
                        )
                    ]),
                    bash.md5sum(
                        copy_done_file + ".md5",
                        None,
                        check=True
                    ),
                    bash.touch(copy_done_file)
                ],
                output_dependency=copy_job_output_dependency + [re.sub(os.path.dirname(self.bioinfo_file), raw_fastq_dir, self.bioinfo_file), copy_done_file],
                name="index.copy_raw." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            )
            jobs.append(copy_job)

        # If copy was already made and successful
        else:
            log.info("Copy of source run folder already done and successful... skipping index.copy_raw." + self.run_id + "." + str(self.lane_number) + " job..." )

        # Then process the copied fastq
        if isinstance(fastq_job, list):
            for job in fastq_job:
                job.input_files.append(copy_done_file)
                jobs.append(job)
        else: 
            jobs.append(fastq_job)

        self.add_copy_job_inputs(jobs)
        return jobs

    def fastq(self):
        """
        *** In the future, may generate the fastq files from the raw CAL files. ***
        When fastq files are already generated by the sequencer, then copy all the fastq files in processing folder.
        The fastq (and metrics) files which match the readsets of the sample sheet are put in the Unaligned folder and renamed accordingly to their record in the sheet.
        All the fastq (and metrics) files which don't match any of the readsets in the sample sheet are put in the Unexpected folder without renaming.
        """
        jobs = []

        if len(self.readsets) == 1:
            log.info("Only one sample on lane : no need to demultiplex...  Skipping fastq step...")

        elif self.is_demultiplexed:
            log.info("Demultiplexing done on the sequencer... Skipping fastq step...")

        else:
            log.info("No demultiplexing done on the sequencer... Processing fastq step...")

            input_fastq_dir = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number), "raw_fastq")

            fastq_multx_outputs, final_fastq_jobs = self.generate_fastqdemultx_outputs() 
            tmp_output_dir = os.path.dirname(fastq_multx_outputs[0])

            fastq_multx_outputs.insert(0, os.path.join(tmp_output_dir, self.flowcell_id + "_L0" + str(self.lane_number) + ".demultiplex.metrics"))
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(tmp_output_dir),
                        run_processing_tools.fastqmultx(
                            os.path.join(self.output_dir, "barcodesheet." + str(self.lane_number) + ".tsv"),
                            self.number_of_mismatches,
                            fastq_multx_outputs,
                            os.path.join(input_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_R1.fastq.gz"),
                            os.path.join(input_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_I1.fastq.gz"),
                            os.path.join(input_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_R2.fastq.gz"),
                            os.path.join(input_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_I2.fastq.gz") if self.is_dual_index else None
                        )
                    ],
                    name="fastq.demultiplex" + self.run_id + "." + str(self.lane_number),
                    samples=self.samples
                )
            )

            if final_fastq_jobs:
                jobs.extend(final_fastq_jobs)

        self.add_copy_job_inputs(jobs)
        return jobs

    def demuxfastq(self):
        """
        *** In the future, may generate the fastq files from the raw CAL files. ***
        When fastq files are already generated by the sequencer, then copy all the fastq files in processing folder.
        The fastq (and metrics) files which match the readsets of the sample sheet are put in the Unaligned folder and renamed accordingly to their record in the sheet.
        All the fastq (and metrics) files which don't match any of the readsets in the sample sheet are put in the Unexpected folder without renaming.
        """
        jobs = []

        if len(self.readset) == 1:
            log.info("Only one sample on lane : no need to demultiplex...  Skipping fastq step...")

        elif self.is_demultiplexed:
            log.info("Demultiplexing done on the sequencer... Skipping fastq step...")

        elif "(no_barcode)" in self.run_dir:
            log.info("No barcode in fastq files... Demultiplexing is undoable... Skipping fastq step...")

        else:
            log.info("No demultiplexing done on the sequencer... Processing fastq step...")

            input_fastq_dir = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number), "raw_fastq")

            demux_fastqs_outputs, final_fastq_jobs = self.generate_demuxfastqs_outputs()
            tmp_output_dir = os.path.dirname(demux_fastqs_outputs[0])

            metrics_file = os.path.join(tmp_output_dir, self.run_id + "." + str(self.lane_number) + ".DemuxFastqs.metrics.txt")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(tmp_output_dir),
                        run_processing_tools.demux_fastqs(
                            os.path.join(self.output_dir, "samplesheet." + str(self.lane_number) + ".csv"),
                            self.number_of_mismatches,
                            self.mask,
                            demux_fastqs_outputs,
                            metrics_file,
                            os.path.join(input_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_1.fq.gz"),
                            os.path.join(input_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz")
                        )
                    ],
                    name="fastq.demultiplex." + self.run_id + "." + str(self.lane_number),
                    samples=self.samples
                )
            )

            if final_fastq_jobs:
                jobs.extend(final_fastq_jobs)

        self.add_copy_job_inputs(jobs)
        return jobs

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
            self.report_inputs['qc'][readset.name] = os.path.join(output_dir, "mpsQC_" + region_name + "_stats.xml")

        self.add_to_report_hash("qc_graphs", jobs)
        self.add_copy_job_inputs(jobs)
        return self.throttle_jobs(jobs)

    def fastqc(self):
        """
        """
        jobs = []
        for readset in self.readsets:
            input_dict = {
                readset.fastq1 : [
                    os.path.join(os.path.dirname(readset.fastq1), "fastqc.R1", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.fastq1))),
                    os.path.join(os.path.dirname(readset.fastq1), "fastqc.R1", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.fastq1)))
                ]
            }
            if readset.run_type == "PAIRED_END":
                input_dict[readset.fastq2] = [
                    os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.fastq2))),
                    os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.fastq2))),
                ]
            for input, outputs in input_dict.items():
                job_suffix = re.search('.*fastqc\.([IR][12])', os.path.dirname(outputs[0])).group(1)
                jobs.append(
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

        self.add_to_report_hash("fastqc", jobs)
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
            self.report_inputs['blast'][readset.name] = os.path.join(output)
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
                name="blast." + readset.name + ".blast." + self.run_id + "." + str(self.lane_number),
                samples = [readset.sample]
            )
            jobs.append(job)

        self.add_to_report_hash("blast", jobs)
        self.add_copy_job_inputs(jobs)
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
        for readset in [readset for readset in self.readsets if readset.bam]:
            job = readset.aligner.get_alignment_job(readset)
            job.samples = [readset.sample]
            jobs.append(job)

        self.add_to_report_hash("align", jobs)
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
            metrics_file = input_file_prefix + "dup.metrics"
            self.report_inputs['mark_dup'][readset.name] = metrics_file

            job = picard.mark_duplicates(
                [input],
                output,
                metrics_file
            )
            job.name = "picard_mark_duplicates." + readset.name + ".dup." + self.run_id + "." + str(self.lane_number)
            job.samples = [readset.sample]
            jobs.append(job)

        self.add_to_report_hash("picard_mark_duplicates", jobs)
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

            if readset.is_rna:
                self.report_inputs['align'][readset.name] = [
                    readset.bam + '.metrics.alignment_summary_metrics',
                    readset.bam + '.metrics.insert_size_metrics',
                    os.path.join(os.path.dirname(readset.bam), readset.sample.name + "." + readset.library + '.rnaseqc.sorted.dup.metrics.tsv'),
                    readset.bam + '.metrics.rRNA.tsv'
                ]
            else:
                self.report_inputs['align'][readset.name] = [
                    readset.bam + '.metrics.alignment_summary_metrics',
                    readset.bam + '.metrics.insert_size_metrics',
                    readset.bam + ".dup.metrics",
                    readset.bam + ".metrics.verifyBamId.tsv",
                    readset.bam + ".metrics.targetCoverage.txt"
                ]

        self.add_to_report_hash("metrics", jobs)
        self.add_copy_job_inputs(jobs)
        return self.throttle_jobs(jobs)

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

            job = concat_jobs(
                current_jobs,
                name="md5." + readset.name + ".md5." + self.run_id + "." + str(self.lane_number),
                samples=[readset.sample]
            )

            jobs.append(job)

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
        Copy the whole processing foler to where they can be serve or loaded into a LIMS
        """
        copy_job = concat_jobs(
            [
                bash.mkdir(
                    os.path.join(
                        config.param("copy", "mgi_runs_root", type="dirpath"),
                        self.seqtype,
                        self.year
                    )
                ),
                bash.cp(
                    os.path.join(self.output_dir, "."),
                    os.path.join(
                        config.param("copy", "mgi_runs_root", type="dirpath"),
                        self.seqtype,
                        self.year,
                        self.instrument + "_" + self.flowcell_id + "_" + self.run_id
                    ),
                    recursive=True
                )
            ],
            input_dependency=self.copy_job_inputs,
            name="copy." + self.run_id + "." + str(self.lane_number),
            samples=self.samples
        )
        return [copy_job]

    def report(self):
        """
        Generate a JSON file reporting the whole pipeline
        """
        jobs = []

        report_dir = os.path.join(self.output_dir, "report")
        sample_report_dir  = os.path.join(report_dir, "sample_json")

        run_validation_inputs = []

        # Add barcodes info to the report_hash
        if self.is_demultiplexed or len(self.readsets) == 1:
            index_per_readset = {}
            for readset in self.readsets:
                index_per_readset[readset.name] = [{
                    'INDEX_NAME': readset.index_name,
                    'INDEX1': None,
                    'INDEX2': None
                }]
            self._index_per_readset = index_per_readset

        self.report_hash["barcodes"] = self.index_per_readset

        general_information_file = os.path.join(self.output_dir, self.run_id + "." + str(self.lane_number) + ".general_information.json")
        with open(general_information_file, 'w') as out_json:
            json.dump(self.report_hash, out_json, indent=4)

        # Build JSON report for each sample
        for readset in self.readsets:
            # Build JSON report for each sample
            output_file = os.path.join(sample_report_dir, readset.name + ".report.json")
            jobs.append(
                concat_jobs([
                    bash.mkdir(sample_report_dir),
                    tools.run_validation_sample_report(
                        readset,
                        self.report_inputs,
                        output_file,
                        general_information_file
                    )],
                    name="sample_report." + readset.name + "." + self.run_id + "." + str(readset.lane)
                )
            )
            run_validation_inputs.append(output_file)

        run_validation_report_json = os.path.join(report_dir, self.run_id + "." + str(self.lane_number) + ".run_validation_report.json")
        # Aggregate sample reports to build the run validation report
        jobs.append(
            concat_jobs([
                bash.mkdir(report_dir),
                tools.run_validation_aggregate_report(
                    general_information_file,
                    run_validation_inputs,
                    run_validation_report_json
                )],
                name="report." + self.run_id + "." + str(self.lane_number)
            )
        )

        # Copy fastqc HTML files into the report folder
        copy_fastqc_jobs = [
            bash.mkdir(os.path.join(self.output_dir, "report", "fastqc"))
        ]
        for readset in self.readsets:
            copy_fastqc_jobs.append(
                concat_jobs([
                    bash.cp(
                        os.path.join(os.path.dirname(readset.fastq1), "fastqc.R1", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.fastq1))),
                        os.path.join(self.output_dir, "report", "fastqc/")
                    ),
                    bash.cp(
                        os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.fastq2))),
                        os.path.join(self.output_dir, "report", "fastqc/")
                )]
            )
        )
        copy_fastqc_jobs[0].output_files.append(os.path.join(self.output_dir, "report", "fastqc"))

        # Copy the MGI summary HTML report into the report folder
        summary_report_html = os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + readset.lane + ".summaryReport.html")
        jobs.append(
            concat_jobs(
                copy_fastqc_jobs + [
                bash.cp(
                    summary_report_html,
                    os.path.join(self.output_dir, "report")
                )],
                name="report.copy." + self.run_id + "." + str(self.lane_number)
            )
        )

        # Create a report zip archive containing all the above...
        zip_output = os.path.join(self.output_dir, "report", self.run_id + "_" + self.flowcell_id + "_L00" + str(self.lane_number) + "_report.zip")
        zip_job =  bash.zip(
            [
                run_validation_report_json,
                os.path.join(self.output_dir, "report", "fastqc"),
                summary_report_html
            ],
            zip_output,
            recursive=False
        )
        zip_job.name = "report.zip." + self.run_id + "." + str(self.lane_number)
        zip_job.samples = self.samples
        jobs.append(zip_job)

        return jobs

    #
    # Utility methods
    #

    def add_copy_job_inputs(self, jobs):
        for job in jobs:
            # we first remove dependencies of the current job, since we will have a dependency on that job
            self.copy_job_inputs = [item for item in self.copy_job_inputs if item not in job.input_files]
            self.copy_job_inputs.extend(job.output_files)

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

    def get_read1cycles(self):
        """
        Parse the BioInfo.csv file for the number of read 1 cycles in the run
        """
        bioinfo_csv = csv.reader(open(self.bioinfo_file, 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Read1 Cycles":
                return row[1]
        _raise(SanitycheckError("Could not get Read 1 Cycles from " + os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv")))

    def get_read2cycles(self):
        """
        Parse the BioInfo.csv file for the number of read 2 cycles in the run
        """
        bioinfo_csv = csv.reader(open(self.bioinfo_file, 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Read2 Cycles":
                return row[1]
        _raise(SanitycheckError("Could not get Read 2 Cycles from " + os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv")))

    def get_index1cycles(self):
        """
        Parse the BioInfo.csv file for the number of index 1 cycles in the run
        """
        bioinfo_csv = csv.reader(open(self.bioinfo_file, 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Barcode":
                return row[1]
        _raise(SanitycheckError("Could not get Index 1 Cycles from " + os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv")))

    def get_index2cycles(self):
        """
        Parse the BioInfo.csv file for the number of index 2 cycles in the run
        """
        bioinfo_csv = csv.reader(open(self.bioinfo_file, 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Dual Barcode":
                if row[1]:
                    return row[1]
                else:
                    return "0"
        _raise(SanitycheckError("Could not get Index 2 Cycles from " + os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv")))

    def get_instrument(self):
        """
        Parse the BioInfo.csv file for the instrument name the run has been running on
        """
        bioinfo_csv = csv.reader(open(self.bioinfo_file, 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Machine ID":
                return row[1]
        _raise(SanitycheckError("Could not find intrument from " + os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv")))

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

    def validate_barcodes(self):
        """
        Validate all index sequences against each other to ensure they aren't in collision according to the chosen
        number of mismatches parameter.
        """
        min_allowed_distance = (2 * self.number_of_mismatches) + 1
        index_lengths = self.get_smallest_index_length()

        validated_indexes = []
        collisions = []

        for readset in self.readsets:
            for current_index in readset.indexes:
                if self.is_dual_index:
                    current_barcode = current_index['INDEX2'][0:index_lengths[0]]+current_index['INDEX1'][0:index_lengths[1]]
                else:
                    current_barcode = current_index['INDEX1'][0:index_lengths[0]]
                for candidate_index in validated_indexes:
                    if distance(current_barcode, candidate_index) < min_allowed_distance:
                        collisions.append("'" + current_barcode + "' and '" + candidate_index + "'")
                validated_indexes.append(current_barcode)

        if len(collisions) > 0:
            _raise(SanitycheckError("Barcode collisions: " + ";".join(collisions)))

    def get_mask(self):
        """
        Returns an fgbio DemuxFastqs friendly mask of the reads cycles.

        The mask is calculated using:
            - first base and last base of index;
            - the index length in the sample sheet;
            - the number of index cycles on the sequencer;
        """
        mask = ""
        index_lengths = self.get_smallest_index_length()
        index_read_count = 0
        nb_total_index_base_used = 0

        index_cycles = [int(self.get_index1cycles())]
        if self.is_dual_index:
            index_cycles.insert(0, int(self.get_index2cycles()))

        mask = self.get_read1cycles() + 'T ' + self.get_read2cycles() + 'T'
        for idx_nb_cycles in index_cycles:
            if idx_nb_cycles >= index_lengths[index_read_count]:
                if index_lengths[index_read_count] == 0 or self.last_index <= nb_total_index_base_used:
                    # Don't use any index bases for this read
                    mask += str(idx_nb_cycles) + 'S'
                else:
                    nb_s_printed = 0

                    # Ns in the beginning of the index read
                    if self.first_index > (nb_total_index_base_used + 1):
                        nb_s_printed = min(idx_nb_cycles, self.first_index - nb_total_index_base_used - 1)
                        if nb_s_printed >= index_lengths[index_read_count]:
                            nb_s_printed = idx_nb_cycles
                        mask += str(nb_s_printed) + 'S'

                    # Calculate the number of index bases
                    nb_index_bases_used = max(index_lengths[index_read_count] - nb_s_printed, 0)
                    nb_index_bases_used = min(self.last_index - nb_total_index_base_used - nb_s_printed, nb_index_bases_used)
                    nb_total_index_base_used += nb_index_bases_used + min(nb_s_printed, index_lengths[index_read_count])
                    if nb_index_bases_used > 0:
                        mask += str(nb_index_bases_used) + 'B'

                    # Ns at the end of the index read
                    remaining_base_count = idx_nb_cycles - nb_index_bases_used - nb_s_printed
                    if remaining_base_count > 0:
                        mask += str(remaining_base_count) + 'S'
            index_read_count += 1

        return mask

    def generate_mgi_lane_sample_sheet(self):
        """
        Create a sample sheet for fgbio DemuxFastqs
        """

        csv_headers = [
            "Sample_ID",
            "Sample_Name",
            "Library_ID",
            "Description",
            "Sample_Barcode"
        ]
        csv_file = os.path.join(
            self.output_dir,
            "samplesheet." + str(self.lane_number) + ".csv"
        )
        writer = csv.DictWriter(
            open(csv_file, 'wb'),
            delimiter=str(','),
            fieldnames=csv_headers
        )

        writer.writeheader()

        # barcode validation
        if re.search("B", self.mask):
            self.validate_barcodes()

        index_lengths = self.get_smallest_index_length()
        for readset in self.readsets:
            for readset_index in readset.indexes:
                # Barcode sequence should only match with the barcode cycles defined in the mask
                # so we adjust thw lenght of the index sequences accordingly for the "Sample_Barcode" field
                if self.is_dual_index:
                    sample_barcode = readset_index['INDEX2'][0:index_lengths[0]] + readset_index['INDEX1'][0:index_lengths[1]]
                else: 
                    sample_barcode = readset_index['INDEX1'][0:index_lengths[0]]
                if self.last_index < len(sample_barcode):
                    sample_barcode = sample_barcode[0:self.last_index]
                if self.first_index > 1:
                    sample_barcode = sample_barcode[self.first_index-1:]

                csv_dict = {
                    "Sample_ID": readset_index['SAMPLESHEET_NAME'],
                    "Sample_Name": readset_index['SAMPLESHEET_NAME'] + '_' + readset_index['INDEX_NAME'],
                    "Library_ID": readset_index['LIBRARY'],
                    "Description": readset.name + '_' + readset.library_type + '_' + readset.library_source,
                    "Sample_Barcode": readset_index['INDEX2_RAW'] + readset_index['INDEX1_RAW']
                }
                writer.writerow(csv_dict)

        self._index_per_readset = index_per_readset

    def generate_fastqdemultx_outputs(self):
        fastq_multx_outputs = []
        final_fastq_jobs = []

        output_dir = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number))

        for readset in self.readsets:
            final_fastq_jobs_per_readset = [
                bash.mkdir(os.path.dirname(readset.fastq1))
            ]

            # If 10X libraries : 4 indexes per sample
            if re.search("SI-", readset.index_name):
                fastqr1a = os.path.join(output_dir, "tmp", readset.name + "_A_R1.fastq")
                fastqr1b = os.path.join(output_dir, "tmp", readset.name + "_B_R1.fastq")
                fastqr1c = os.path.join(output_dir, "tmp", readset.name + "_C_R1.fastq")
                fastqr1d = os.path.join(output_dir, "tmp", readset.name + "_D_R1.fastq")
                fastq_multx_outputs.extend([
                     fastqr1a,
                     fastqr1b,
                     fastqr1c,
                     fastqr1d,
                     undet_fastqr1 if undet_fastqr1 else None
                ])

                final_fastq_jobs_per_readset.append(
                    pipe_jobs([
                        bash.cat(
                            list(filter(None, [
                                fastqr1a,
                                fastqr1b,
                                fastqr1c,
                                fastqr1d,
                                undet_fastqr1 if undet_fastqr1 else None
                            ])),
                            None
                        ),
                        bash.gzip(
                            None,
                            readset.fastq1
                        )
                    ])
                )

                # Add the fastq of first index
                fastqi1a = os.path.join(output_dir, "tmp", readset.name + "_A_I1.fastq")
                fastqi1b = os.path.join(output_dir, "tmp", readset.name + "_B_I1.fastq")
                fastqi1c = os.path.join(output_dir, "tmp", readset.name + "_C_I1.fastq")
                fastqi1d = os.path.join(output_dir, "tmp", readset.name + "_D_I1.fastq")
                fastq_multx_outputs.extend([
                    fastqi1a,
                    fastqi1b,
                    fastqi1c,
                    fastqi1d,
                    undet_fastqi1 if undet_fastqi1 else None
                ])

                final_fastq_jobs_per_readset.append(
                    pipe_jobs(
                        [
                            bash.cat(
                                list(filter(None, [
                                    fastqi1a,
                                    fastqi1b,
                                    fastqi1c,
                                    fastqi1d,
                                    undet_fastqi1 if undet_fastqi1 else None
                                ])),
                                None
                           ),
                           bash.gzip(
                                None,
                                readset.index_fastq1
                           )
                        ]
                    )
                )

                # For paired-end sequencing, do not forget the fastq of the reverse reads
                if readset.run_type == "PAIRED_END" :
                    fastqr2a = os.path.join(output_dir, "tmp", readset.name + "_A_R2.fastq")
                    fastqr2b = os.path.join(output_dir, "tmp", readset.name + "_B_R2.fastq")
                    fastqr2c = os.path.join(output_dir, "tmp", readset.name + "_C_R2.fastq")
                    fastqr2d = os.path.join(output_dir, "tmp", readset.name + "_D_R2.fastq")
                    fastq_multx_outputs.extend([
                        fastqr2a,
                        fastqr2b,
                        fastqr2c,
                        fastqr2d,
                        undet_fastqr2 if undet_fastqr2 else None
                    ])

                    final_fastq_jobs_per_readset.append(
                        pipe_jobs([
                            bash.cat(
                                list(filter(None, [
                                    fastqr2a,
                                    fastqr2b,
                                    fastqr2c,
                                    fastqr2d,
                                    undet_fastqr2 if undet_fastqr2 else None
                                ])),
                                None
                            ),
                            bash.gzip(
                                None,
                                readset.fastq2
                            )
                        ])
                    )

                    # For dual index demultiplexing, do not forget the fastq of the second index
                    if self.is_dual_index:
                        fastqi2a = os.path.join(output_dir, "tmp", readset.name + "_A_I2.fastq")
                        fastqi2b = os.path.join(output_dir, "tmp", readset.name + "_B_I2.fastq")
                        fastqi2c = os.path.join(output_dir, "tmp", readset.name + "_C_I2.fastq")
                        fastqi2d = os.path.join(output_dir, "tmp", readset.name + "_D_I2.fastq")
                        fastq_multx_outputs.extend([
                            fastqi2a,
                            fastqi2b,
                            fastqi2c,
                            fastqi2d,
                            undet_fastqi2 if undet_fastqi2 else None
                        ])

                        final_fastq_jobs_per_readset.append(
                            pipe_jobs([
                                bash.cat(
                                    list(filter(None, [
                                        fastqi2a,
                                        fastqi2b,
                                        fastqi2c,
                                        fastqi2d,
                                        undet_fastqi2 if undet_fastqi2 else None
                                    ])),
                                    None
                                ),
                                bash.gzip(
                                    None,
                                    readset.index_fastq2
                                )
                            ])
                        )

            # not a 10X library : 1 index per sample
            else:
                fastqr1 = os.path.join(output_dir, "tmp", readset.name + "_R1.fastq")
                fastq_multx_outputs.extend([
                    fastqr1,
                    undet_fastqr1 if undet_fastqr1 else None
                ])
                final_fastq_jobs_per_readset.append(
                    pipe_jobs([
                        bash.cat(
                            list(filter(None, [
                                fastqr1,
                                undet_fastqr1 if undet_fastqr1 else None
                            ])),
                            None
                        ),
                        bash.gzip(
                            None,
                            readset.fastq1
                        )
                    ])
                )
                fastqi1 = os.path.join(output_dir, "tmp", readset.name + "_I1.fastq")
                fastq_multx_outputs.extend([
                    fastqi1,
                    undet_fastqi1 if undet_fastqi1 else None
                ])
                final_fastq_jobs_per_readset.append(
                    pipe_jobs([
                        bash.cat(
                            list(filter(None, [
                                fastqi1,
                                undet_fastqi1 if undet_fastqi1 else None
                            ])),
                            None
                        ),
                        bash.gzip(
                            None,
                            readset.index_fastq1
                        )
                    ])
                )
                # For paired-end sequencing, do not forget the fastq of the reverse reads
                if readset.run_type == "PAIRED_END":
                    fastqr2 = os.path.join(output_dir, "tmp", readset.name + "_R2.fastq")
                    fastq_multx_outputs.extend([
                        fastqr2,
                        undet_fastqr2 if undet_fastqr2 else None
                    ])
                    final_fastq_jobs_per_readset.append(
                        pipe_jobs([
                            bash.cat(
                                list(filter(None, [
                                    fastqr2,
                                    undet_fastqr2 if undet_fastqr2 else None
                                ])),
                                None
                            ),
                            bash.gzip(
                                None,
                                readset.fastq2
                            )
                        ])
                    )
                    # For dual index multiplexing, do not forget the fastq of the second index
                    if self.is_dual_index:
                        fastqi2 = os.path.join(output_dir, "tmp", readset.name + "_I2.fastq")
                        fastq_multx_outputs.extend([
                            fastqi2,
                            undet_fastqi2 if undet_fastqi2 else None
                        ])
                        final_fastq_jobs_per_readset.append(
                            pipe_jobs([
                                bash.cat(
                                    list(filter(None, [
                                        fastqi2,
                                        undet_fastqi2 if undet_fastqi2 else None
                                    ])),
                                    None
                                ),
                                bash.gzip(
                                    None,
                                    readset.index_fastq2
                                )
                            ])
                        )

            if len(final_fastq_jobs_per_readset) != 0:
                final_fastq_jobs.append(
                    concat_jobs(
                        final_fastq_jobs_per_readset,
                        name="fastq.rename."+readset.name,
                        samples=[readset.sample]
                    )
                )

        return fastq_multx_outputs, final_fastq_jobs

    def generate_demuxfastqs_outputs(self):
        demuxfastqs_outputs = []
        final_fastq_jobs = []
        count = 0

        output_dir = os.path.join(self.output_dir, "Unaligned.test4." + str(self.lane_number))

        for readset in self.readsets:
            readset_r1_outputs = []
            readset_r2_outputs = []
            readset_r1_rename_jobs = []
            readset_r2_extract_jobs = []
            readset_i1_extract_jobs = []
            readset_i2_extract_jobs = []

#            final_fastq_jobs_per_readset = [
#                bash.mkdir(os.path.dirname(readset.fastq1))
#            ]
            readset_r1_rename_jobs.append(
                bash.mkdir(os.path.dirname(readset.fastq1))
            )

            for index in readset.indexes:
                readset_r1_outputs.append(
                    os.path.join(output_dir, "tmp", index['SAMPLESHEET_NAME']+"-"+index['SAMPLESHEET_NAME']+'_'+index['INDEX_NAME']+"-"+index['INDEX2_RAW']+index['INDEX1_RAW']+"_R1.fastq.gz")
                )
                if readset.run_type == "PAIRED_END":
                    readset_r2_outputs.append(
                        os.path.join(output_dir, "tmp", index['SAMPLESHEET_NAME']+"-"+index['SAMPLESHEET_NAME']+'_'+index['INDEX_NAME']+"-"+index['INDEX2_RAW']+index['INDEX1_RAW']+"_R2.fastq.gz")
                    )

            # If True, then merge the 'Undetermined' reads
            if self.merge_undetermined:
                readset_r1_outputs.append(
                    os.path.join(output_dir, "tmp", "unmatched_R1.fastq.gz")
                )
                if readset.run_type == "PAIRED_END":
                    readset_r2_outputs.append(
                        os.path.join(output_dir, "tmp", "unmatched_R2.fastq.gz")
                    )

            if len(readset_r1_outputs) > 1:
                readset_rename_jobs.append(
                    pipe_jobs(
                        [
                            bash.cat(
                                readset_r1_outputs,
                                None,
                                zip=True
                            ),
                            bash.pigz(
                                None,
                                readset.fastq1
                            )
                        ]
                    )
                )
            else:
                readset_r1_rename_jobs.append(
                    bash.mv(
                        readset_r1_outputs[0],
                        readset.fastq1
                    )
                )

            if readset.run_type == "PAIRED_END":
                demuxfastqs_outputs.extend(readset_r2_outputs)
                raw_readset_fastq2 = os.path.join(
                    os.path.dirname(readset_r2_outputs[0]),
                    re.sub("_001.fastq.gz", "_001.raw.fastq.gz", os.path.basename(readset.fastq2))
                )
                # if more than one barcode used for the readset, aggregate all R2 fastq into one file for the current readset
                if len(readset_r2_outputs) > 1:
                    readset_rename_jobs.append(
                        pipe_jobs(
                            [
                                bash.cat(
                                    readset_r2_outputs,
                                    None,
                                    zip=True
                                ),
                                bash.pigz(
                                    None,
                                    raw_readset_fastq2
                                )
                            ]
                        )
                    )
                # if only one barcode, then just rename the file
                else:
                    readset_rename_jobs.append(
                        bash.mv(
                            readset_r2_outputs[0],
                            raw_readset_fastq2
                        )
                    )

                final_fastq_jobs_per_readset.append(
                    bash.cat(
                        readset_r1_outputs,
                        None,
                        zip=True
                    ),
                    bash.gzip(
                        None,
                        readset.fastq1
                    )
                )
                if readset.run_type == "PAIRED_END":
                    final_fastq_jobs_per_readset.append(
                        bash.cat(
                            readset_r2_outputs,
                            None,
                            zip=True
                        ),
                        bash.cut(
                            None,
                            None,
                            options="-c 1-" + self.get_read2cycles()
                        ),
                        bash.pigz(
                            None,
                            readset.fastq2
                        )
                    ],
                    name="fastq.extract_R2_fastq." + readset.name + "." + self.run_id + "." + str(self.lane_number),
                    samples=self.samples
                )
            )
            if self.is_dual_index:
                # I1 fastq
                final_fastq_jobs.append(
                    pipe_jobs(
                        [
                            bash.cat(
                                raw_readset_fastq2,
                                None,
                                zip=True
                            ),
                            bash.paste(
                                None,
                                None,
                                options="-d '\\t' - - - -"
                            ),
                            bash.awk(
                                None,
                                None,
                                "-F'\\t' '{print $1 \"\\n\" substr($2,"+str(int(self.get_read2cycles())+int(self.get_index2cycles())+1)+","+self.get_index1cycles()+") \"\\n\" $3 \"\\n\" substr($4,"+str(int(self.get_read2cycles())+int(self.get_index2cycles())+1)+","+self.get_index1cycles()+") }'"
                            ),
                            bash.pigz(
                                None,
                                readset.index_fastq1
                            )
                        ],
                        name="fastq.extract_I1_fastq." + readset.name + "." + self.run_id + "." + str(self.lane_number)
                    )

            if len(readset_demuxfastqs_outputs) > 1:
                final_fastq_jobs_per_readset.append(
                    bash.cat(
                        list(filter(None, readset_demuxfastqs_outputs)),
                        None,
                        zip=True
                    ),
                    bash.gzip(
                        None,
                        readset.fastq1
                    )
                )
            else:
                final_fastq_jobs_per_readset.append(
                    
                )
                # I2 fastq
                final_fastq_jobs.append(
                    pipe_jobs(
                        [
                            bash.cat(
                                list(filter(None, [
                                    fastqr2a,
                                    fastqr2b,
                                    fastqr2c,
                                    fastqr2d,
                                    undet_fastqr2 if undet_fastqr2 else None
                                ])),
                                None
                            ),
                            bash.pigz(
                                None,
                                readset.fastq2
                            )
                        ])
                    )

            # not a 10X library : 1 index per sample
            else:
                fastqr1 = os.path.join(output_dir, "tmp", readset.name + "_R1.fastq")
                demuxfastqs_outputs.extend([
                    fastqr1,
                    undet_fastqr1 if undet_fastqr1 else None
                ])
                final_fastq_jobs_per_readset.append(
                    pipe_jobs([
                        bash.cat(
                            list(filter(None, [
                                fastqr1,
                                undet_fastqr1 if undet_fastqr1 else None
                            ])),
                            None
                        ),
                        bash.gzip(
                            None,
                            readset.fastq1
                        )
                    ])
                )
                # For paired-end sequencing, do not forget the fastq of the reverse reads
                if readset.run_type == "PAIRED_END":
                    fastqr2 = os.path.join(output_dir, "tmp", readset.name + "_R2.fastq")
                    demuxfastqs_outputs.extend([
                        fastqr2,
                        undet_fastqr2 if undet_fastqr2 else None
                    ])
                    final_fastq_jobs_per_readset.append(
                        pipe_jobs([
                            bash.cat(
                                list(filter(None, [
                                    fastqr2,
                                    undet_fastqr2 if undet_fastqr2 else None
                                ])),
                                None
                            ),
                            bash.pigz(
                                None,
                                readset.fastq2
                            )
                        ])
                    )

            if len(final_fastq_jobs_per_readset) != 0:
                final_fastq_jobs.append(
                    concat_jobs(
                        final_fastq_jobs_per_readset,
                        name="fastq.rename."+readset.name,
                        samples=[readset.sample]
                    )
                )
        return demuxfastqs_outputs, final_fastq_jobs

    def get_smallest_index_length(self):
        """
        Returns a list (for each index read of the run) of the minimum between the number of index cycle on the
        sequencer and all the index lengths.
        """
        run_index_lengths = []

        all_indexes = []
        for readset in self.readsets:
            all_indexes += readset.index

        if self.is_dual_index:
            min_sample_index_length = min(len(index['INDEX2']) for index in all_indexes)
            run_index_lengths.append(min(min_sample_index_length, int(self.get_index2cycles())))

        min_sample_index_length = min(len(index['INDEX1']) for index in all_indexes)
        run_index_lengths.append(min(min_sample_index_length, int(self.get_index1cycles())))

        return run_index_lengths

    def load_readsets(self):
        """
        Parse the sample sheet and return a list of readsets.
        """

        return parse_mgi_raw_readset_files(
            self.readset_file,
            self.bioinfo_file,
            "PAIRED_END" if self.is_paired_end else "SINGLE_END",
            self.seqtype,
            self.flowcell_id,
            int(self.lane_number),
            self.get_read1cycles(),
            self.get_read2cycles(),
            self.output_dir
        )

    def submit_jobs(self):
        super(MGIRunProcessing, self).submit_jobs()

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
            self.index,
#            self.fastq,
            self.demuxfastq,
            self.qc_graphs,
            self.fastqc,
            self.blast,
            self.align,
            self.picard_mark_duplicates,
            self.metrics,
        #    self.md5,
            self.copy,
            self.report
        ]

def distance(str1, str2):
    """
    Returns the hamming distance. http://code.activestate.com/recipes/499304-hamming-distance/#c2
    """
    return sum(itertools.imap(str.__ne__, str1, str2))

if __name__ == '__main__':

    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        MGIRunProcessing()

