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
        self.argparser.add_argument("-x", help="First index base to use for demultiplexing (inclusive). The index from the sample sheet will be adjusted according to that value.", type=int, required=False, dest="first_index")
        self.argparser.add_argument("-y", help="Last index base to use for demultiplexing (inclusive)", type=int, required=False, dest="last_index")
        self.argparser.add_argument("-m", help="Number of index mistmaches allowed for demultiplexing (default 1). Barcode collisions are always checked.", type=int, required=False, dest="number_of_mismatches")

        super(MGIRunProcessing, self).__init__(protocol)

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = self.load_readsets()
            if not self.is_demultiplexed:
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

    def is_dual_index(self):
        if not hasattr(self, "_is_dual_index"):
            if self.get_index2cycles():
                self._is_dual_index = True
            else:
                self._is_dual_index = False
        return self._is_dual_index

    @property
    def is_demultiplexed(self):
        if not hasattr(self, "_is_demultiplexed"):
            if all(readset.is_mgi_index for readset in self.readsets):
                self._is_demultiplexed = True
            elif all(not readset.is_mgi_index for readset in self.readsets):
                self._is_demultiplexed = False
            else:
                _raise(SanitycheckError("Error: bad index settings in lane... both non-MGI and MGI adapters were detected in the readset file" + self.readset_file))
        return self._is_demultiplexed

    @property
    def mask(self):
        if not hasattr(self, "_mask"):
            if self.is_dual_index:
                self._mask = self.get_read1cycles() + "T " + self.get_read2cycles() + "T" + self.get_index2cycles() + "B" + self.get_index1cycles() + "B"
            else:
                self._mask = self.get_read1cycles() + "T " + self.get_read2cycles() + "T" + self.get_index1cycles() + "B"
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
        """
        The instrument id from the run folder
        """
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
            #for readset in self.readsets:
            #    log.error(readset.run)
            runs = set([readset.run for readset in self.readsets])
            #log.error(runs)
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
    def flowcell_id(self):
        """
        The flow cell ID from the run folder
        """
        if not hasattr(self, "_flowcell_id"):
            self._flowcell_id = os.path.basename(self.run_dir.rstrip('/'))
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
        copy_job_output_dependency = []
        fastq_job_input_dependency = []

        # If demultiplexing were done on the sequencer i.e. MGI adpaters used...
        if self.is_demultiplexed:
            fastq_jobs = []

            # ...then do the necessary moves and renamings, from the raw_fastq folder to the Unaligned folders
            for readset in self.readsets:
                output_dir = os.path.join(self.output_dir, "Unaligned." + readset.lane, 'Project_' + readset.project, 'Sample_' + readset.name)

                copy_job_output_dependency.extend([
                    readset.fastq1,
                    re.sub("gz", "fqStat.txt", readset.fastq1),
                    re.sub("_R1.fastq.gz", ".report.html", readset.fastq1)
                ])
                fastq_job_input_dependency.extend([
                    os.path.join(self.run_dir, "L0" + str(self.lane_number), readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_1.fq.gz"),
                    os.path.join(self.run_dir, "L0" + str(self.lane_number), readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_1.fq.fqStat.txt"),
                    os.path.join(self.run_dir, "L0" + str(self.lane_number), readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + ".report.html")
                ])
                fastq_jobs.append(
                    concat_jobs([
                        bash.mkdir(output_dir),
                        bash.mv(
                            os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_1.fq.gz"),
                            readset.fastq1
                        ),
                        bash.ln(
                            readset.fastq1,
                            os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_1.fq.gz")
                        ),
                        bash.mv(
                            os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_1.fq.fqStat.txt"),
                            re.sub("gz", "fqStat.txt", readset.fastq1)
                        ),
                        bash.ln(
                            re.sub("gz", "fqStat.txt", readset.fastq1),
                            os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_1.fq.fqStat.txt")
                        ),
                        bash.mv(
                            os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + ".report.html"),
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
                        readset.fastq2,
                        re.sub("gz", "fqStat.txt", readset.fastq2)
                    ])
                    fastq_job_input_dependency.extend([
                        os.path.join(self.run_dir, "L0" + str(self.lane_number), readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_2.fq.gz"),
                        os.path.join(self.run_dir, "L0" + str(self.lane_number), readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_2.fq.fqStat.txt")
                    ])
                    fastq_jobs.append(
                        concat_jobs([
                            bash.mv(
                                os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_2.fq.gz"),
                                readset.fastq2,
                            ),
                            bash.ln(
                                readset.fastq2,
                                os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_2.fq.gz")
                            ),
                            bash.mv(
                                os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index_number + "_2.fq.fqStat.txt"),
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

        # If no demultiplexing were done on the sequencer i.e. no MGI adpater used... 
        else:
            # ...then extract read2, index1 and index2 from R2 fastq, and do some renamings before demultiplexing ('fastq' step)

            # R1 fastq
            copy_job_output_dependency.append(os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_R1.fastq.gz"))
            fastq_job_input_dependency.append(os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_1.fq.gz"))
            fastq_job = [concat_jobs(
                [
                    bash.mv(
                        os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_1.fq.gz"),
                        os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_R1.fastq.gz")
                    ),
                    bash.ln(
                        os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_R1.fastq.gz"),
                        os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_read_1.fq.gz")
                    )
                ],
                input_dependency=[os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_1.fq.gz")],
                output_dependency=[os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_R1.fastq.gz")],
                name="index.rename_R1_fastq." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            )]
            if self.is_paired_end:
                # R2 fastq
                copy_job_output_dependency.append(os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_R2.fastq.gz"))
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
                                os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_R2.fastq.gz"),
                            )
                        ],
                        input_dependency=[os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz")],
                        name="index.extract_R2_fastq." + self.run_id + "." + str(self.lane_number),
                        samples=self.samples
                    )
                )
            if self.is_dual_index:
                # I1 fastq
                copy_job_output_dependency.append(os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_I1.fastq.gz"))
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
                                os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_I1.fastq.gz"),
                            )
                        ],
                        input_dependency=[os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz")],
                        name="index.extract_I1_fastq." + self.run_id + "." + str(self.lane_number),
                        samples=self.samples
                    )
                )
                # I2 fastq
                copy_job_output_dependency.append(os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_I2.fastq.gz"))
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
                                os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_I2.fastq.gz"),
                            )
                        ],
                        input_dependency=[os.path.join(self.run_dir, "L0" + str(self.lane_number), self.flowcell_id + "_L0" + str(self.lane_number) + "_read_2.fq.gz")],
                        name="index.extract_I2_fastq." + self.run_id + "." + str(self.lane_number),
                        samples=self.samples
                    )
                )
            else:
                # I1 fastq
                copy_job_output_dependency.append(os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_I1.fastq.gz"))
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
                                os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_I1.fastq.gz"),
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
                            list(set(fastq_job_input_dependency)),
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
                        ),
                    ]),
                    bash.md5sum(
                        copy_done_file + ".md5",
                        None,
                        check=True
                    ),
                    bash.touch(copy_done_file)
                ],
                output_dependency=copy_job_output_dependency + [copy_done_file],
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

        if self.is_demultiplexed:
            log.info("Demultiplexing done on the sequencer... Skipping fastq step...")

        else:
            log.info("No demultiplexing done on the sequencer... Processing fastq step...")

            input_fastq_dir = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number), "raw_fastq")

            fastq_multx_outputs, final_fastq_jobs = self.generate_fastqdemultx_outputs() 
            tmp_output_dir = os.path.dirname(fastq_multx_outputs[0])

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

        if self.is_demultiplexed:
            log.info("Demultiplexing done on the sequencer... Skipping fastq step...")

        else:
            log.info("No demultiplexing done on the sequencer... Processing fastq step...")

            input_fastq_dir = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number), "raw_fastq")

            demux_fastqs_outputs, final_fastq_jobs = self.generate_demuxfastqs_outputs()
            tmp_output_dir = os.path.dirname(demux_fastqs_outputs[0])

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(tmp_output_dir),
                        run_processing_tools.demux_fastqs(
                            os.path.join(self.output_dir, "samplesheet." + str(self.lane_number) + ".csv"),
                            self.number_of_mismatches,
                            self.mask,
                            demux_fastqs_outputs,
                            os.path.join(input_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_R1.fastq.gz"),
                            os.path.join(input_fastq_dir, self.flowcell_id + "_L0" + str(self.lane_number) + "_R2.fastq.gz")
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
#            if not self.is_demultiplexed:
#                input_dict[readset.index_fastq1] = [
#                    os.path.join(os.path.dirname(readset.index_fastq1), "fastqc.I1", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.index_fastq1))),
#                    os.path.join(os.path.dirname(readset.index_fastq1), "fastqc.I1", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.index_fastq1)))
#                ]
#                if readset.index_type == "DUALINDEX":
#                    input_dict[readset.index_fastq2] = [
#                        os.path.join(os.path.dirname(readset.index_fastq2), "fastqc.I2", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.index_fastq2))),
#                        os.path.join(os.path.dirname(readset.index_fastq2), "fastqc.I2", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.index_fastq2))),
#                    ]
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

#            self.report_inputs['index'].append(readset.bam + '.metrics.alignment_summary_metrics')

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
        if self.is_demultiplexed:
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
            zip_output
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
        bioinfo_csv = csv.reader(open(os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv"), 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Read1 Cycles":
                return row[1]
        _raise(SanitycheckError("Could not get Read 1 Cycles from " + os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv")))

    def get_read2cycles(self):
        """
        Parse the BioInfo.csv file for the number of read 2 cycles in the run
        """
        bioinfo_csv = csv.reader(open(os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv"), 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Read2 Cycles":
                return row[1]
        _raise(SanitycheckError("Could not get Read 2 Cycles from " + os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv")))

    def get_index1cycles(self):
        """
        Parse the BioInfo.csv file for the number of index 1 cycles in the run
        """
        bioinfo_csv = csv.reader(open(os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv"), 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Barcode":
                return row[1]
        _raise(SanitycheckError("Could not get Index 1 Cycles from " + os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv")))

    def get_index2cycles(self):
        """
        Parse the BioInfo.csv file for the number of index 2 cycles in the run
        """
        bioinfo_csv = csv.reader(open(os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv"), 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Dual Barcode":
                return row[1]
        _raise(SanitycheckError("Could not get Index 2 Cycles from " + os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv")))

    def get_instrument(self):
        """
        Parse the BioInfo.csv file for the instrument name the run has been running on
        """
        bioinfo_csv = csv.reader(open(os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv"), 'rb'))
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

    def generate_mgi_lane_barcode_sheet(self):
        """
        Create a barcode sheet to use with the fastq-multix software
        from the Clarity data.

        Only the samples of the chosen lane will be in the file.
        The sample indexes are trimmed according to the mask used.
        """

        csv_file = os.path.join(
            self.output_dir,
            "barcodesheet." + str(self.lane_number) + ".tsv"
        )
        writer = csv.writer(open(csv_file, 'wb'), dialect="excel-tab")

        index_per_readset = {}
        for readset in self.readsets:
            readset_indexes = self.get_index(readset)
            index_per_readset[readset.name] = readset_indexes

            for readset_index in readset_indexes:
                if readset.index_type == "DUALINDEX":
                    index_to_use = readset_index['INDEX1']+"-"+ readset_index['INDEX2']
                else:
                    index_to_use = readset_index['INDEX1']
                csv_data = [
                    readset_index['BCL2FASTQ_NAME'],
                    index_to_use,
                    readset.index_name
                ]
                writer.writerow(csv_data)

        self._index_per_readset = index_per_readset

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

        count = 0
        index_per_readset = {}
        for readset in self.readsets:
            count += 1
            readset_indexes = self.get_index(readset)
            index_per_readset[readset.name] = readset_indexes

            for readset_index in readset_indexes:

                csv_dict = {
                    "Sample_ID": "Sample_" + readset_index['BCL2FASTQ_NAME'],
                    "Sample_Name": readset_index['BCL2FASTQ_NAME'],
                    "Library_ID": readset_index['LIBRARY'],
                    "Description": readset.name + ' - ' + readset.library_type + ' - ' + readset.library_source,
                    "Sample_Barcode": readset_index['INDEX1'] + "-" +  readset_index['INDEX2']
                }
                writer.writerow(csv_dict)

        self._index_per_readset = index_per_readset

    def generate_fastqdemultx_outputs(self):
        fastq_multx_outputs = []
        final_fastq_jobs = []
        count = 0

        output_dir = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number))

        # If True, then merge the 'Undetermined' reads
        if self.merge_undetermined:
            undet_fastqr1 = os.path.join(output_dir, "tmp", "unmatched_R1.fastq")
            undet_fastqi1 = os.path.join(output_dir, "tmp", "unmatched_I1.fastq")
            undet_fastqr2 = os.path.join(output_dir, "tmp", "unmatched_R2.fastq") if readset.run_type == "PAIRED_END" else None
            undet_fastqi2 = os.path.join(output_dir, "tmp", "unmatched_I2.fastq") if readset.index_type == "DUALINDEX" else None
        else:
            undet_fastqr1 = None
            undet_fastqi1 = None
            undet_fastqr2 = None
            undet_fastqi2 = None

        for readset in self.readsets:
            final_fastq_jobs_per_readset = [
                bash.mkdir(os.path.dirname(readset.fastq1))
            ]

            fastq_output_dir = os.path.join(output_dir, 'Project_' + readset.project, 'Sample_' + readset.name)
            # If 10X libraries : 4 indexes per sample
            if re.search("tenX", readset.library_type):
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
                    if readset.index_type == "DUALINDEX" :
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
                    if readset.index_type == "DUALINDEX" :
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

        output_dir = os.path.join(self.output_dir, "Unaligned.test." + str(self.lane_number))

        # If True, then merge the 'Undetermined' reads
        if self.merge_undetermined:
            undet_fastqr1 = os.path.join(output_dir, "tmp", "unmatched_R1.fastq")
            undet_fastqi1 = os.path.join(output_dir, "tmp", "unmatched_I1.fastq")
            undet_fastqr2 = os.path.join(output_dir, "tmp", "unmatched_R2.fastq") if readset.run_type == "PAIRED_END" else None
            undet_fastqi2 = os.path.join(output_dir, "tmp", "unmatched_I2.fastq") if readset.index_type == "DUALINDEX" else None
        else:
            undet_fastqr1 = None
            undet_fastqi1 = None
            undet_fastqr2 = None
            undet_fastqi2 = None

        for readset in self.readsets:
            final_fastq_jobs_per_readset = [
                bash.mkdir(os.path.dirname(readset.fastq1))
            ]

            fastq_output_dir = os.path.join(output_dir, 'Project_' + readset.project, 'Sample_' + readset.name)
            # If 10X libraries : 4 indexes per sample
            if re.search("tenX", readset.library_type):
                fastqr1a = os.path.join(output_dir, "tmp", readset.name + "_A_R1.fastq")
                fastqr1b = os.path.join(output_dir, "tmp", readset.name + "_B_R1.fastq")
                fastqr1c = os.path.join(output_dir, "tmp", readset.name + "_C_R1.fastq")
                fastqr1d = os.path.join(output_dir, "tmp", readset.name + "_D_R1.fastq")
                demuxfastqs_outputs.extend([
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

#                # Add the fastq of first index
#                fastqi1a = os.path.join(output_dir, "tmp", readset.name + "_A_I1.fastq")
#                fastqi1b = os.path.join(output_dir, "tmp", readset.name + "_B_I1.fastq")
#                fastqi1c = os.path.join(output_dir, "tmp", readset.name + "_C_I1.fastq")
#                fastqi1d = os.path.join(output_dir, "tmp", readset.name + "_D_I1.fastq")
#                fastq_multx_outputs.extend([
#                    fastqi1a,
#                    fastqi1b,
#                    fastqi1c,
#                    fastqi1d,
#                    undet_fastqi1 if undet_fastqi1 else None
#                ])
#
#                final_fastq_jobs_per_readset.append(
#                    pipe_jobs(
#                        [
#                            bash.cat(
#                                list(filter(None, [
#                                    fastqi1a,
#                                    fastqi1b,
#                                    fastqi1c,
#                                    fastqi1d,
#                                    undet_fastqi1 if undet_fastqi1 else None
#                                ])),
#                                None
#                           ),
#                           bash.gzip(
#                                None,
#                                readset.index_fastq1
#                           )
#                        ]
#                    )
#                )

                # For paired-end sequencing, do not forget the fastq of the reverse reads
                if readset.run_type == "PAIRED_END" :
                    fastqr2a = os.path.join(output_dir, "tmp", readset.name + "_A_R2.fastq")
                    fastqr2b = os.path.join(output_dir, "tmp", readset.name + "_B_R2.fastq")
                    fastqr2c = os.path.join(output_dir, "tmp", readset.name + "_C_R2.fastq")
                    fastqr2d = os.path.join(output_dir, "tmp", readset.name + "_D_R2.fastq")
                    demuxfastqs_outputs.extend([
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

#                    # For dual index demultiplexing, do not forget the fastq of the second index
#                    if readset.index_type == "DUALINDEX" :
#                        fastqi2a = os.path.join(output_dir, "tmp", readset.name + "_A_I2.fastq")
#                        fastqi2b = os.path.join(output_dir, "tmp", readset.name + "_B_I2.fastq")
#                        fastqi2c = os.path.join(output_dir, "tmp", readset.name + "_C_I2.fastq")
#                        fastqi2d = os.path.join(output_dir, "tmp", readset.name + "_D_I2.fastq")
#                        fastq_multx_outputs.extend([
#                            fastqi2a,
#                            fastqi2b,
#                            fastqi2c,
#                            fastqi2d,
#                            undet_fastqi2 if undet_fastqi2 else None
#                        ])
#
#                        final_fastq_jobs_per_readset.append(
#                            pipe_jobs([
#                                bash.cat(
#                                    list(filter(None, [
#                                        fastqi2a,
#                                        fastqi2b,
#                                        fastqi2c,
#                                        fastqi2d,
#                                        undet_fastqi2 if undet_fastqi2 else None
#                                    ])),
#                                    None
#                                ),
#                                bash.gzip(
#                                    None,
#                                    readset.index_fastq2
#                                )
#                            ])
#                        )

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
#                fastqi1 = os.path.join(output_dir, "tmp", readset.name + "_I1.fastq")
#                fastq_multx_outputs.extend([
#                    fastqi1,
#                    undet_fastqi1 if undet_fastqi1 else None
#                ])
#                final_fastq_jobs_per_readset.append(
#                    pipe_jobs([
#                        bash.cat(
#                            list(filter(None, [
#                                fastqi1,
#                                undet_fastqi1 if undet_fastqi1 else None
#                            ])),
#                            None
#                        ),
#                        bash.gzip(
#                            None,
#                            readset.index_fastq1
#                        )
#                    ])
#                )
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
                            bash.gzip(
                                None,
                                readset.fastq2
                            )
                        ])
                    )
#                    # For dual index multiplexing, do not forget the fastq of the second index
#                    if readset.index_type == "DUALINDEX" :
#                        fastqi2 = os.path.join(output_dir, "tmp", readset.name + "_I2.fastq")
#                        demuxfastqs_outputs.extend([
#                            fastqi2,
#                            undet_fastqi2 if undet_fastqi2 else None
#                        ])
#                        final_fastq_jobs_per_readset.append(
#                            pipe_jobs([
#                                bash.cat(
#                                    list(filter(None, [
#                                        fastqi2,
#                                        undet_fastqi2 if undet_fastqi2 else None
#                                    ])),
#                                    None
#                                ),
#                                bash.gzip(
#                                    None,
#                                    readset.index_fastq2
#                                )
#                            ])
#                        )

            if len(final_fastq_jobs_per_readset) != 0:
                final_fastq_jobs.append(
                    concat_jobs(
                        final_fastq_jobs_per_readset,
                        name="fastq.rename."+readset.name,
                        samples=[readset.sample]
                    )
                )

        return demuxfastqs_outputs, final_fastq_jobs

    def get_index(self, readset):
        """
        """

        indexes = []
        idx_seq_array = []
        index1 = ""
        index2 = ""
        index1seq = ""
        index2seq = ""

        index_file = config.param('DEFAULT', 'index_settings_file', type='filepath', required=False)
        if not (index_file and os.path.isfile(index_file)):
            index_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_settings_format.txt')

        index_seq_pattern = "grep '%s,' %s | head -n1 | awk -F',' '{print $%d}'"

        if readset.library_type == "tenX_sc_RNA_v2" or readset.library_type == "tenX_DNA_v2" or readset.library_type == "FeatureBarcode" or readset.library_type == "tenX_sc_ATAC_v1":
            index1 = readset.index_name

            for idx, char in enumerate(['_A', '_B', '_C', '_D']):
                index1seq = subprocess.check_output(index_seq_pattern % (index1, index_file, idx+2), shell=True).strip()
                [actual_index1seq, actual_index2seq, adapteri7, adapteri5] = self.sub_get_index(readset, index1seq, index2seq)
                indexes.append({
                    'BCL2FASTQ_NAME': readset.name + char,
                    'LIBRARY': readset.library,
                    'PROJECT': readset.project,
                    'INDEX_NAME': readset.index_name,
                    'INDEX1': actual_index1seq,
                    'INDEX2': actual_index2seq,
                    'ADAPTERi7' : adapteri7,
                    'ADAPTERi5' : adapteri5
                })
                (actual_index1seq != "") and idx_seq_array.append(actual_index1seq)
                (actual_index2seq != "") and idx_seq_array.append(actual_index2seq)

        elif readset.library_type == "tenX_sc_RNA_v1" :
            index2 = readset.index_name


            for idx, char in enumerate(['_A', '_B', '_C', '_D']):
                index2seq = subprocess.check_output(index_seq_pattern % (index2, index_file, idx+2), shell=True).strip()
                [actual_index1seq, actual_index2seq, adapteri7, adapteri5] = self.sub_get_index(readset, index1seq, index2seq)
                indexes.append({
                    'BCL2FASTQ_NAME': readset.name + char,
                    'LIBRARY': readset.library,
                    'PROJECT': readset.project,
                    'INDEX_NAME': readset.index_name,
                    'INDEX1': actual_index1seq,
                    'INDEX2': actual_index2seq,
                    'ADAPTERi7' : adapteri7,
                    'ADAPTERi5' : adapteri5
                })
                (actual_index1seq != "") and idx_seq_array.append(actual_index1seq)
                (actual_index2seq != "") and idx_seq_array.append(actual_index2seq)

        else:
            if re.search("-", readset.index_name):
                index1 = readset.index_name.split("-")[0]
                index2 = readset.index_name.split("-")[1]
                index1seq = subprocess.check_output(index_seq_pattern % (index1, index_file, 2), shell=True).strip()
                index2seq = subprocess.check_output(index_seq_pattern % (index2, index_file, 2), shell=True).strip()

            else:
                index1 = readset.index_name

                index1seq = subprocess.check_output(index_seq_pattern % (index1, index_file, 2), shell=True).strip()

            [actual_index1seq, actual_index2seq, adapteri7, adapteri5] = self.sub_get_index(readset, index1seq, index2seq)
            indexes.append({
                'BCL2FASTQ_NAME': readset.name,
                'LIBRARY': readset.library,
                'PROJECT': readset.project,
                'INDEX_NAME': readset.index_name,
                'INDEX1': actual_index1seq,
                'INDEX2': actual_index2seq,
                'ADAPTERi7' : adapteri7,
                'ADAPTERi5' : adapteri5
            })
            (actual_index1seq != "") and idx_seq_array.append(actual_index1seq)
            (actual_index2seq != "") and idx_seq_array.append(actual_index2seq)

        return indexes

    def sub_get_index(self, readset, index1seq, index2seq):
        """
        """
        index_file = config.param('DEFAULT', 'index_settings_file', type='filepath', required=False)
        if not (index_file and os.path.isfile(index_file)):
            index_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_settings_format.txt')

        index1cycles = int(self.get_index1cycles())
        index2cycles = int(self.get_index2cycles())

        actual_index1seq='';
        actual_index2seq='';

        seqtype = self.seqtype
        primer_seq_pattern = "grep -A8 '%s' %s | grep '_IDX_' | awk -F':' '{print $2}' | tr -d \"35\'\- \" | awk -F',' '{print $%d}'"

        index1_primer = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index 1") % (seqtype, index_file, 1), shell=True).strip()
        index1_primeroffset = int(subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index 1") % (seqtype, index_file, 2), shell=True).strip())
        index2_primer = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index 2") % (seqtype, index_file, 1), shell=True).strip()
        index2_primeroffset = int(subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index 2") % (seqtype, index_file, 2), shell=True).strip())

        main_seq_pattern = "grep -A4 '^%s:' %s | grep -E \"^3'|^5'\" | head -n 1 | sed \"s/5'//g\"  | sed \"s/3'//g\" | tr -d \" '\-\" | grep %s"
        actual_seq_pattern = "echo %s | awk -F\"%s\" '{print $%d}' | sed \"%s\" | cut -c %s"

        present = subprocess.check_output(main_seq_pattern % (readset.library_type, index_file, "-c " + index1_primer), shell=True).strip()
        if present == "1" :

            main_seq = subprocess.check_output(main_seq_pattern % (readset.library_type, index_file, index1_primer), shell=True).strip()

            if seqtype == "dnbseqg400" :
                actual_index1seq = subprocess.check_output(actual_seq_pattern % (main_seq, index1_primer, 2, "s/\[i7\]/"+index1seq+"/g", str(index1_primeroffset+1)+"-"+str(index1_primeroffset+index1cycles)+" | tr 'ATGC' 'TACG' | rev"), shell=True).strip()
            else:
                actual_index1seq = subprocess.check_output(actual_seq_pattern % (main_seq, index1_primer, 2, "s/\[i7\]/"+index1seq+"/g", str(index1_primeroffset+1)+"-"+str(index1_primeroffset+index1cycles)), shell=True).strip()

            if readset.index_type == "DUALINDEX" :
                main_seq = subprocess.check_output(main_seq_pattern.replace("| head -n 1 |", "|") % (readset.library_type, index_file, index2_primer), shell=True).strip()

                if seqtype == "hiseqx" or seqtype == "hiseq4000" or seqtype == "iSeq" or seqtype == "dnbseqg400" :
                    actual_index2seq = subprocess.check_output(actual_seq_pattern.replace("| cut -c", "| rev | cut -c") % (main_seq, index2_primer, 1, "s/\[i5c\]/$(echo "+index2seq+" | tr 'ATGC' 'TACG' )/g", str(index2_primeroffset+1)+"-"+str(index2_primeroffset+index2cycles)), shell=True).strip()
                else :
                    actual_index2seq = subprocess.check_output(actual_seq_pattern % (main_seq, index2_primer, 2, "s/\[i5\]/"+index2seq+"/g", str(index2_primeroffset+1)+"-"+str(index2_primeroffset+index2cycles)), shell=True).strip()

        else :
            indexn1_primer = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index N1").replace("_DIGIT_", "1"), shell=True).strip()
            indexn1_primeroffset = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index N1").replace("_DIGIT_", "2"), shell=True).strip()
            indexn2_primer = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index N2").replace("_DIGIT_", "1"), shell=True).strip()
            indexn2_primeroffset = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index N2").replace("_DIGIT_", "2"), shell=True).strip()

            main_seq = subprocess.check_output(main_seq_pattern % (readset.library_type, index_file, indexn1_primer), shell=True).strip()

            if seqtype == "dnbseqg400" :
                actual_index1seq = subprocess.check_output(actual_seq_pattern % (main_seq, indexn1_primer, 2, "s/\[i7\]/"+index1seq+"/g", str(indexn1_primeroffset+1)+"-"+str(indexn1_primeroffset+index1cycles))+" | tr 'ATGC' 'TACG' | rev", shell=True).strip()
            else:
                actual_index1seq = subprocess.check_output(actual_seq_pattern % (main_seq, indexn1_primer, 2, "s/\[i7\]/"+index1seq+"/g", str(indexn1_primeroffset+1)+"-"+str(indexn1_primeroffset+index1cycles)), shell=True).strip()

            if readset.index_type == "DUALINDEX" :
                main_seq = subprocess.check_output(main_seq_pattern.replace("| head -n 1 |", "|") % (readset.library_type, index_file, indexn2_primer), shell=True).strip()

                if seqtype == "hiseqx" or seqtype == "hiseq4000" or seqtype == "iSeq" or seqtype == "dnbseqg400" :
                    actual_index2seq = subprocess.check_output(actual_seq_pattern.replace("| cut -c", "| rev | cut -c") % (main_seq, indexn2_primer, 1, "s/\[i5c\]/$(echo "+index2seq+" | tr 'ATGC' 'TACG' )/g", str(indexn2_primeroffset+1)+"-"+str(indexn2_primeroffset+index2cycles)), shell=True).strip()
                else :
                    actual_index2seq = subprocess.check_output(actual_seq_pattern % (main_seq, indexn2_primer, 2, "s/\[i5\]/"+index2seq+"/g", str(indexn2_primeroffset+1)+"-"+str(indexn2_primeroffset+index2cycles)), shell=True).strip()

        main_seq = subprocess.check_output(main_seq_pattern % (readset.library_type, index_file, "\[i7\]"), shell=True).strip()
        adapteri7 = subprocess.check_output("echo \"%s\" | awk -F\'\\\[i7\\\]\' '{print $1}' | awk -F\'\\\]\' '{print $NF}'" % (main_seq), shell=True).strip()

        main_seq = subprocess.check_output(main_seq_pattern.replace("| head -n 1 |", "|") % (readset.library_type, index_file, "\[i5c\]"), shell=True).strip()
        adapteri5 = subprocess.check_output("echo \"%s\" | awk -F\'\\\[i5c\\\]\' '{print $2}' | awk -F\'\\\[\' '{print $1}' | rev" % (main_seq), shell=True).strip()

        return [actual_index1seq, actual_index2seq, adapteri7, adapteri5]


    def load_readsets(self):
        """
        Parse the sample sheet and return a list of readsets.
        """

        return parse_mgi_raw_readset_files(
            self.readset_file,
            self.run_dir,
            "PAIRED_END" if self.is_paired_end else "SINGLE_END",
            self.flowcell_id,
            int(self.lane_number),
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
            self.fastq,
#            self.demuxfastq,
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

if __name__ == '__main__':

    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        MGIRunProcessing()
