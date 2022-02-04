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
#from __future__ import print_function, division, unicode_literals, absolute_import
import argparse
import logging
import os
import re
import sys
import shutil
import xml.etree.ElementTree as Xml
import math
import csv
from collections import OrderedDict, namedtuple
import json

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs

from bfx.readset import parse_mgi_raw_readset_files
from bfx import run_processing_tools
from bfx import bash_cmd as bash

from pipelines import run_processing_common as common

log = logging.getLogger(__name__)

class MGIRunProcessing(common.RunProcessing):
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
        self.copy_job_inputs = {}
        self.argparser.add_argument("-f", "--flag", help="T7 flag files directory (mandatory for T7 runs)", required='type' in sys.argv, dest="raw_flag_dir")
        self.argparser.add_argument("--raw-fastq-prefix", help="Prefix used to search for the raw fastq from the sequencer. Default <FLOWCELL_ID>_<RUN_ID>", required=False, dest="raw_fastq_prefix")
        self.argparser.add_argument("--splitbarcode-demux", help="demultiplexing done while basecalling with MGI splitBarcode", action="store_true", required=False, dest="splitbarcode_demux")
        self.argparser.add_argument("-t", "--type", help = "MGI sequencing technology", choices = ["g400", "t7"], default="g400")

        args = sys.argv[1:]
        t7 = False
        flag = False
        for i, arg in enumerate(args):
            if arg in ['-t', '--type'] and args[i+1] == 't7':
                t7 = True
            if arg in ['-f', '--flag']:
                flag = True
        if flag and not t7:
            _raise(SanitycheckError("-f/--flag option can only be set when -t/--type=t7."))

        super(MGIRunProcessing, self).__init__(protocol)

    @property
    def run_dir(self):
        if self.args.run_dir:
            return self.args.run_dir
        else:
            _raise(SanitycheckError("Error: missing '-d/--run_dir' option!"))

    @property
    def raw_flag_dir(self):
        if self.args.raw_flag_dir:
            return self.args.raw_flag_dir
        else:
            _raise(SanitycheckError("Error: missing '-f/--flag' option!"))

    @property
    def flowcell_id(self):
        """
        The flow cell ID from the run folder or from parameter
        """
        if not hasattr(self, "_flowcell_id"):
            if self.args.flowcell_id:
                self._flowcell_id = self.args.flowcell_id
            else:
                rundir_basename = os.path.basename(self.run_dir.rstrip('/'))
                if "_" in rundir_basename:
                    [self._flowcell_id, junk_food] = rundir_basename.split("_")
                else:
                    self._flowcell_id = os.path.basename(self.run_dir.rstrip('/'))
                    if "(" in self._flowcell_id:
                         self._flowcell_id = self._flowcell_id.split("(")[0]
        return self._flowcell_id

    @property
    def run_id(self):
        """
        The RUN ID from the run folder or from parameter
        """
        if not hasattr(self, "_run_id"):
            if self.args.run_id:
                self._run_id = self.args.run_id
            else:
                rundir_basename = os.path.basename(self.run_dir.rstrip('/'))
                if "_" in rundir_basename:
                    [junk_food, self._run_id] = rundir_basename.split("_")
                else:
                    _raise(SanitycheckError("Error: Run ID could not be parsed from the RUN folder : " + self.run_dir))
        return self._run_id

    @property
    def raw_fastq_prefix(self):
        """
        """
        if not hasattr(self, "_raw_fastq_prefix"):
            if self.args.raw_fastq_prefix:
                self._raw_fastq_prefix = self.args.raw_fastq_prefix
            else:
                self._raw_fastq_prefix = self.flowcell_id + "_" + self.run_id
        return self._raw_fastq_prefix

    @property
    def bioinfo_files(self):
        if not hasattr(self, "_bioinfo_files"):
            self._bioinfo_files = {}
            for lane in self.lanes:
                if os.path.exists(os.path.join(self.output_dir, "Unaligned." + lane, "raw_fastq", "BioInfo.csv")):
                    self._bioinfo_files[lane] = os.path.join(self.output_dir, "Unaligned." + lane, "raw_fastq", "BioInfo.csv")
                else:
                    bioinfo_file = os.path.join(self.run_dir, "L0" + lane, "BioInfo.csv")
                    if lane in self._bioinfo_files and self._bioinfo_files[lane] != bioinfo_file:
                        _raise(SanitycheckError("More than one Bioinfo.csv found for lane '" + lane + "' : " + self._bioinfo_files[lane] + ", " + bioinfo_file))
                    else:
                        self._bioinfo_files[lane] = bioinfo_file
        return self._bioinfo_files

    @property
    def no_index_fastq(self):
        if not hasattr(self, "_no_index_fastq"):
            if self.args.splitbarcode_demux:
                self._no_index_fastq = True 
            else:
                self._no_index_fastq = False
        return self._no_index_fastq

    @property
    def bioinfo_hash(self):
        if not hasattr(self, "_bioinfo_hash"):
            self._bioinfo_hash = {}
            for lane in self.lanes:
                bioinfo_csv = csv.reader(open(self.bioinfo_files[lane], 'r'))
                self._bioinfo_hash[lane] = dict(bioinfo_csv)
        return self._bioinfo_hash

    @property
    def json_flag_files(self):
        if not hasattr(self, "_json_flag_files"):
            self._json_flag_files = {}
            for lane in self.lanes:
                self._json_flag_files[lane] = self.generate_mgi_t7_flag_file(lane)
        return self._json_flag_files

    @property
    def json_flag_hash(self):
        if not hasattr(self, "_json_flag_hash"):
            self._json_flag_hash = {}
            for lane in self.lanes:
                with open(self.json_flag_files[lane], "r") as jff:
                    json_flag_content = json.load(jff)
                keys = [item[0] for item in json_flag_content['experimentInfoVec']]
                vals = [item[1] for item in json_flag_content['experimentInfoVec']]
                self._json_flag_hash[lane] = dict(zip(keys, vals))
        return self._json_flag_hash

    def basecall(self):
        """
        Use write_fastq software from MGI to perfor the base calling.
        Takes the raw .cal files from the sequencer and produces fastq files.
        Demultiplexing while doing the basecalling is still under testings...
        """

        jobs = []

        for lane in self.lanes:
            lane_jobs = []
            jobs_to_throttle = []

            input = self.readset_file

            unaligned_dir = os.path.join(self.output_dir, "Unaligned." + lane)
            basecall_dir = os.path.join(unaligned_dir, "basecall")

            lane_config_file = os.path.join(unaligned_dir, self.run_id + "." + lane + ".settings.config")

            # If demultiplexing is perform while basecalling
            if self.args.splitbarcode_demux:
                basecall_outputs, postprocessing_jobs = self.generate_basecall_outputs(lane)
                basecall_outputs.extend(
                    [
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + ".summaryReport.html"),
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + ".heatmapReport.html"),
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, "summaryTable.csv")
                    ]
                )

                lane_jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(basecall_dir),
                            run_processing_tools.mgi_t7_basecall(
                                input,
                                self.flowcell_id,
                                basecall_outputs,
                                basecall_dir,
                                self.json_flag_files[lane],
                                lane_config_file
                            )
                        ],
                        name="basecall." + self.run_id + "." + lane,
                        samples=self.samples[lane]
                    )
                )

                if postprocessing_jobs:
                    jobs_to_throttle.extend(postprocessing_jobs)

                self.report_inputs[lane]['basecall'] = os.path.join(basecall_dir, self.run_id, "L0" + lane)

            else:
                basecall_outputs = [
                    os.path.join(basecall_dir, self.run_id, "L0" + lane),
                    os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_read_1.fq.gz"),
                    os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_read_2.fq.gz"),
                    os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + ".summaryReport.html"),
                    os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + ".heatmapReport.html"),
                    os.path.join(basecall_dir, self.run_id, "L0" + lane, "summaryTable.csv")
                ]

                raw_fastq_dir = os.path.join(unaligned_dir, "raw_fastq")
                raw_fastq_outputs = [
                    os.path.join(raw_fastq_dir, self.raw_fastq_prefix +  "_L0" + lane + "_read_1.fq.gz"),
                    os.path.join(raw_fastq_dir, self.raw_fastq_prefix +  "_L0" + lane + "_read_2.fq.gz"),
                    os.path.join(raw_fastq_dir, self.raw_fastq_prefix +  "_L0" + lane + ".summaryReport.html"),
                    os.path.join(raw_fastq_dir, self.raw_fastq_prefix +  "_L0" + lane + ".heatmapReport.html"),
                    os.path.join(raw_fastq_dir, "summaryTable.csv")
                ]

                lane_basecall_job = concat_jobs(
                    [
                        bash.mkdir(basecall_dir),
                        run_processing_tools.mgi_t7_basecall(
                            input,
                            self.run_dir,
                            self.flowcell_id,
                            basecall_outputs,
                            basecall_dir,
                            self.json_flag_files[lane],
                            lane_config_file
                        ),
                        bash.ln(
                            os.path.join(basecall_dir, self.run_id, "L0" + lane),
                            raw_fastq_dir
                        )
                    ],
                    name="basecall." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
                lane_basecall_job.output_files.extend(raw_fastq_outputs)
                lane_jobs.append(lane_basecall_job)

            lane_jobs.extend(jobs_to_throttle)

            self.add_to_report_hash("basecall", lane, lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)

        return jobs


    def index(self):
        """
        First copy all the files of the lane from the sequencer deposit folder to the processing folder,
        into "raw_fastq".
        Then, if demultiplexing already done on sequencer, formerly the case when MGI adapters were used, then
            rename the fastq files and move them from "raw_fastq" to "Unaligned.LANE" folder,
            --> pipeline will SKIP demultiplexing (i.e. next step, fastq).
        Else, (i.e. demultiplexing still remains to be done) nothing remains to do here :
            everything is in place for the demultiplexing to happen in the next step
        *TO DO* - in both cases, retrieve index stats from the sequencer output files to build a proper index report
        """

        jobs = []

        for lane in self.lanes:
            if int(self.index1cycles[lane]) + int(self.index2cycles[lane]) == 0 and len(self.readsets[lane]) > 1:
                 err_msg = "LANE SETTING ERROR :\n"
                 err_msg += "Unable to demultiplex " + str(len(self.readsets[lane])) + " samples : No barcode in fastq files...\n(in "
                 err_msg += self.run_dir + ")"
                 _raise(SanitycheckError(err_msg))

            lane_jobs = []

            # First we copy all the lane folder into the raw_fastq folder
            unaligned_dir = os.path.join(self.output_dir, "Unaligned." + lane)
            raw_fastq_dir = os.path.join(unaligned_dir, "raw_fastq")
            copy_done_file = os.path.join(raw_fastq_dir, "copy_done.Success")

            raw_name_prefix = self.raw_fastq_prefix +  "_L0" + lane
            copy_job_output_dependency = [
                os.path.join(raw_fastq_dir, raw_name_prefix + ".summaryReport.html"),
                os.path.join(raw_fastq_dir, raw_name_prefix + ".heatmapReport.html"),
                os.path.join(raw_fastq_dir, "summaryTable.csv")
            ]
            rename_job = None

            # If demultiplexing was done while bawecalling
            if self.args.splitbarcode_demux:
                rename_jobs = []
                copy_job_output_dependency.extend([
                    os.path.join(raw_fastq_dir, "BarcodeStat.txt"),
                    os.path.join(raw_fastq_dir, "SequenceStat.txt")
                ])

                # ...then do the renamings from the raw_fastq folder to the Unaligned folder
                for readset in self.readsets[lane]:
                    output_dir = os.path.join(self.output_dir, "Unaligned." + lane, 'Project_' + readset.project_id, 'Sample_' + readset.name)
                    m = re.search("\D+(?P<index>\d+)", readset.index_name)
                    if m:
                        index_number = m.group('index').lstrip("0")
                    else:
                        _raise(SanitycheckError("MGI Index error : unable to parse barcode # in : " + readset.index_name))

                    copy_job_output_dependency.extend([
                        os.path.join(raw_fastq_dir, raw_name_prefix + "_" + index_number + "_1.fq.gz"),
                        os.path.join(raw_fastq_dir, raw_name_prefix + "_" + index_number + "_1.fq.fqStat.txt"),
                        os.path.join(raw_fastq_dir, raw_name_prefix + "_" + index_number + ".report.html")
                    ])
                    rename_jobs.append(
                        concat_jobs([
                            bash.mkdir(output_dir),
                            bash.cp(
                                os.path.join(raw_fastq_dir, raw_name_prefix + "_" + index_number + "_1.fq.gz"),
                                readset.fastq1
                            ),
                            bash.cp(
                                os.path.join(raw_fastq_dir, raw_name_prefix + "_" + index_number + "_1.fq.fqStat.txt"),
                                re.sub("gz", "fqStat.txt", readset.fastq1)
                            ),
                            bash.cp(
                                os.path.join(raw_fastq_dir, raw_name_prefix + "_" + index_number + ".report.html"),
                                re.sub("_R1_001.fastq.gz", ".report.html", readset.fastq1)
                            )
                        ])
                    )

                    if readset.run_type == "PAIRED_END":
                        copy_job_output_dependency.extend([
                            os.path.join(raw_fastq_dir, raw_name_prefix + "_" + index_number + "_2.fq.gz"),
                            os.path.join(raw_fastq_dir, raw_name_prefix + "_" + index_number + "_2.fq.fqStat.txt")
                        ])
                        rename_jobs.append(
                            concat_jobs([
                                bash.cp(
                                    os.path.join(raw_fastq_dir, raw_name_prefix + "_" + index_number + "_2.fq.gz"),
                                    readset.fastq2,
                                ),
                                bash.cp(
                                    os.path.join(raw_fastq_dir, raw_name_prefix + "_" + index_number + "_2.fq.fqStat.txt"),
                                    re.sub("gz", "fqStat.txt", readset.fastq2),
                                )
                            ])
                        )

                rename_job = concat_jobs(rename_jobs)
                rename_job.output_files.append(copy_done_file)
                rename_job.name = "index.rename." + self.run_id + "." + lane
                rename_job.samples = self.samples[lane]

            # If no demultiplexing was done on the sequencer i.e. no MGI adpater used...
            # Nothing to do, fgbio DemuxFastqs will used the raw fastq files in the fastq step,
            # But set the copy jobs output_dependency
            else:
                copy_job_output_dependency.append(os.path.join(raw_fastq_dir, raw_name_prefix + "_read_1.fq.gz"))
                copy_job_output_dependency.append(os.path.join(raw_fastq_dir, raw_name_prefix + "_read_2.fq.gz"))
                copy_job_output_dependency.append(os.path.join(raw_fastq_dir, raw_name_prefix + "_read_1.fq.fqStat.txt"))
                copy_job_output_dependency.append(os.path.join(raw_fastq_dir, raw_name_prefix + "_read_2.fq.fqStat.txt"))
                copy_job_output_dependency.append(os.path.join(raw_fastq_dir, raw_name_prefix + "_read.report.html"))

                if (len(self.readsets[lane]) == 1) or (int(self._index1cycles[lane]) + int(self._index2cycles[lane]) == 0):
                    readset = self.readsets[lane][0]
                    rename_job = concat_jobs(
                        [
                            bash.mkdir(os.path.dirname(readset.fastq1)),
                            bash.cp(
                                os.path.join(raw_fastq_dir, raw_name_prefix + "_read_1.fq.gz"),
                                readset.fastq1,
                            ),
                            bash.cp(
                                os.path.join(raw_fastq_dir, raw_name_prefix + "_read_1.fq.fqStat.txt"),
                                re.sub("gz", "fqStat.txt", readset.fastq1),
                            ),
                            bash.cp(
                                os.path.join(raw_fastq_dir, raw_name_prefix + "_read.report.html"),
                                re.sub("_R1_001.fastq.gz", "_read.report.html", readset.fastq1),
                            )
                        ]
                    )
                    if readset.run_type == "PAIRED_END":
                        rename_job = concat_jobs(
                            [
                                rename_job,
                                bash.cp(
                                    os.path.join(raw_fastq_dir, raw_name_prefix + "_read_2.fq.gz"),
                                    readset.fastq2,
                                ),
                                bash.cp(
                                    os.path.join(raw_fastq_dir, raw_name_prefix + "_read_2.fq.fqStat.txt"),
                                 re.sub("gz", "fqStat.txt", readset.fastq2),
                                )
                            ]
                        )
                    rename_job.output_files.append(copy_done_file)
                    rename_job.name = "index.rename." + self.run_id + "." + lane
                    rename_job.samples = self.samples[lane]

            # Job submission :
            # First we copy all the lane folder into the raw_fastq folder (unless already done...)
            if not os.path.exists(copy_done_file) or self.force_jobs:
                copy_job = concat_jobs(
                    [
                        bash.mkdir(raw_fastq_dir),
                        bash.cp(
                            os.path.join(self.run_dir, "L0" + lane, "."),
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
                copy_job = concat_jobs(
                    [
                        copy_job,
                        pipe_jobs([
                            bash.md5sum(
                                self.bioinfo_files[lane],
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
                    output_dependency=copy_job_output_dependency + [re.sub(os.path.dirname(self.bioinfo_files[lane]), raw_fastq_dir, self.bioinfo_files[lane]), copy_done_file],
                    name="index.copy_raw." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
                lane_jobs.append(copy_job)

            # If copy was already made and successful
            else:
                log.info("Copy of source run folder already done and successful... skipping \"index.copy_raw." + self.run_id + "." + lane + "\" job..." )

            # Then process the copied fastq
            if rename_job:
                lane_jobs.append(rename_job)

            self.add_to_report_hash("index", lane, lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)

        return jobs

    def fastq(self):
        """
        *** In the future, may generate the fastq files from the raw CAL files. ***
        Perform demultplexing of the reads with fgbio DemuxFastqs
        """

        ini_section = None
        if self.args.type == 'g400':
            ini_section = 'fastq_g400'
        elif self.args.type == 't7':
            ini_section = 'fastq_t7'
        if not ini_section:
            _raise(SanitycheckError("Could not determine which section to use for fastq step from given protocol " + self.protocol))
        jobs = []

        if self.args.splitbarcode_demux:
            log.info("Demultiplexing done during the basecalling... Skipping fastq step...")

        else:
            log.info("No demultiplexing done yet... Processing fastq step...")

            for lane in self.lanes:
                if int(self._index1cycles[lane]) + int(self._index2cycles[lane]) == 0:
                    log.info("No barcode cycles in the lane... Skipping fastq step for lane " + lane + "...")

                else:
                    lane_jobs = []
                    jobs_to_throttle = []

                    input_fastq_dir = os.path.join(self.output_dir, "Unaligned." + lane, "raw_fastq")
                    raw_name_prefix = self.raw_fastq_prefix +  "_L0" + lane
                    input1 = os.path.join(input_fastq_dir, raw_name_prefix + "_read_1.fq.gz")
                    input2 = os.path.join(input_fastq_dir, raw_name_prefix + "_read_2.fq.gz")

                    demuxfastqs_outputs, postprocessing_jobs = self.generate_demuxfastqs_outputs(lane)

                    tmp_output_dir = os.path.dirname(demuxfastqs_outputs[0])
                    tmp_metrics_file = os.path.join(tmp_output_dir, self.run_id + "." + lane + ".DemuxFastqs.metrics.txt")
                    metrics_file = os.path.join(self.output_dir, "Unaligned." + lane, self.run_id + "." + lane + ".DemuxFastqs.metrics.txt")
                    demuxfastqs_outputs.append(metrics_file)

                    if ini_section == 'fastq_g400':
                        demultiplex_job = run_processing_tools.demux_fastqs(
                            os.path.join(self.output_dir, "samplesheet." + lane + ".csv"),
                            self.number_of_mismatches,
                            self.mask[lane],
                            demuxfastqs_outputs,
                            tmp_metrics_file,
                            input1,
                            input2,
                            ini_section
                        )
                    elif ini_section == 'fastq_t7':
                        demultiplex_job = run_processing_tools.demux_fastqs_by_chunk(
                            os.path.join(self.output_dir, "samplesheet." + lane + ".csv"),
                            self.number_of_mismatches,
                            self.mask[lane],
                            demuxfastqs_outputs,
                            tmp_metrics_file,
                            input1,
                            input2,
                            ini_section
                        )
                    lane_jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(
                                    tmp_output_dir,
                                    remove=True
                                ),
                                demultiplex_job,
                                bash.cp(
                                    tmp_metrics_file,
                                    metrics_file
                                )
                            ],
                            name=ini_section + ".demultiplex." + self.run_id + "." + lane,
                            samples=self.samples[lane],
                            input_dependency=demultiplex_job.input_files,
                            output_dependency=demuxfastqs_outputs
                        )
                    )

                    if postprocessing_jobs:
                        jobs_to_throttle.extend(postprocessing_jobs)

                    self.report_inputs[lane]['fastq'] = metrics_file

                    if self.args.type == 't7':
                        lane_jobs.extend(jobs_to_throttle)
                    else:
                        lane_jobs.extend(self.throttle_jobs(jobs_to_throttle))

                    self.add_to_report_hash("fastq", lane, lane_jobs)
                    self.add_copy_job_inputs(lane_jobs, lane)
                    jobs.extend(lane_jobs)

        return jobs

    #
    # Utility methods
    #

    def get_read1cycles(self, lane):
        """
        Parse the BioInfo.csv file for the number of read 1 cycles in the run
        """
        if self.args.type == "g400":
            return self.bioinfo_hash[lane]["Read1 Cycles"]

        elif self.args.type == "t7":
            return self.json_flag_hash[lane]['Read1']

        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_read2cycles(self, lane):
        """
        Parse the BioInfo.csv file for the number of read 2 cycles in the run
        """
        if self.args.type == "g400":
            return self.bioinfo_hash[lane]["Read2 Cycles"]

        elif self.args.type == "t7":
            return self.json_flag_hash[lane]['Read2']
        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_index1cycles(self, lane):
        """
        Parse the BioInfo.csv file for the number of index 1 cycles in the run
        """
        if self.args.type == "g400":
            if self.bioinfo_hash[lane]["Barcode"] == '':
                return "0"
            else:
                return self.bioinfo_hash[lane]["Barcode"]
            _raise(SanitycheckError("Could not get Index 1 Cycles from " + self.bioinfo_files[lane]))

        elif self.args.type == "t7":
            if self.json_flag_hash[lane]["Barcode"] == '':
                return "0"
            else:
                return self.json_flag_hash[lane]['Barcode']
        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_index2cycles(self, lane):
        """
        Parse the BioInfo.csv file for the number of index 2 cycles in the run
        """
        if self.args.type == "g400":
            if self.bioinfo_hash[lane]["Dual Barcode"] == '':
                return "0"
            else:
                return self.bioinfo_hash[lane]["Dual Barcode"]

        elif self.args.type == "t7":
            if self.json_flag_hash[lane]["Dual Barcode"] == '':
                return "0"
            else:
                return self.json_flag_hash[lane]["Dual Barcode"]
        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_instrument(self, lane):
        """
        Parse the BioInfo.csv file for the instrument name the run has been running on
        """
        if self.args.type == "g400":
            return self.bioinfo_hash[lane]["Machine ID"]

        elif self.args.type == "t7":
            return self.json_flag_hash[lane]['Machine ID']
        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_flowcell_position(self, lane):
        """
        Parse the BioInfo.csv file for flowcell position of the run
        """
        if self.args.type == "g400":
            return self.bioinfo_hash[lane]["Flowcell Pos"]

        elif self.args.type == "t7":
            return self.json_flag_hash[lane]["Flow Cell Pos"]

        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_sequencer_run_id(self, lane):
        """
        Parse the BioInfo.csv file for the ID of the run given by the sequencer
        """
        if self.args.type == "g400":
            # dnb_id format looks like : 10074MG01B_Lane4
            # where run_id is : 10074MG01B
            return self.bioinfo_hash[lane]["DNB ID"].split("_")[0]

        elif self.args.type == "t7":
            return self.json_flag_hash[lane]['Flow Cell ID'].split("_")[0]

        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_run_number(self, lane):
        """
        Parse the BioInfo.csv file for the ID of the run to extract the run counter
        """
        if self.args.type == "g400":
            # dnb_id format looks like : 10074MG01B_Lane4
            # where run_id is : 10074MG01B
            # and run counter is : 10074
            run_number = self.bioinfo_hash[lane]["DNB ID"].split("_")[0][:-5]
            # sometimes, format is misleading : 1074
            # so we correct it to : 10074
            while len(run_number) < 5:
                run_number[:1] + "0" +  run_number[1:]
            return run_number

        elif self.args.type == "t7":
            return self.json_flag_hash[lane]['Flow Cell ID'].split("_")[0][-5:]

        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def validate_barcodes(self, lane):
        """
        Validate all index sequences against each other to ensure they aren't in collision according to the chosen
        number of mismatches parameter.
        """
        min_allowed_distance = (2 * self.number_of_mismatches) + 1
        index_lengths = self.get_smallest_index_length(lane)

        validated_indexes = []
        collisions = []

        for readset in self.readsets[lane]:
            for current_index in readset.indexes:
                if self.is_dual_index[lane]:
                    current_barcode = current_index['INDEX2'][0:index_lengths[0]]+current_index['INDEX1'][0:index_lengths[1]]
                else:
                    current_barcode = current_index['INDEX1'][0:index_lengths[0]]
                for candidate_index in validated_indexes:
                    if distance(current_barcode, candidate_index) < min_allowed_distance:
                        collisions.append("'" + current_barcode + "' and '" + candidate_index + "'")
                validated_indexes.append(current_barcode)

        if len(collisions) > 0:
            _raise(SanitycheckError("Barcode collisions: " + ";".join(collisions)))

    def get_mask(self, lane):
        """
        Returns an fgbio DemuxFastqs friendly mask of the reads cycles.

        The mask is calculated using:
            - first base and last base of index;
            - the index length in the sample sheet;
            - the number of index cycles on the sequencer;
        """
        mask = ""
        index_lengths = self.get_smallest_index_length(lane)
        index_read_count = 0
        nb_total_index_base_used = 0

        index_cycles = [int(self.index1cycles[lane])]
        if self.is_dual_index[lane]:
            index_cycles.insert(0, int(self.index2cycles[lane]))

        mask = self.read1cycles[lane] + 'T ' + self.read2cycles[lane] + 'T'
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

    def generate_mgi_t7_flag_file(self, lane):
        """
        Copy the flag file output from the sequencer and if necessary (i.e. t7 fast-mode) insert the barcode information
        This Flag file is used to call splitBarcode
        """

        json_flag_file = os.path.join(self.output_dir, "flag." + lane + ".json")

        raw_flag_file = ""
        for filename in os.listdir(self.raw_flag_dir):
            if re.match(self.run_id + "_" + lane + "_.+json", filename):
                raw_flag_file = os.path.join(self.raw_flag_dir, filename)
                break
        else:
            _raise(SanitycheckError("Could not find any proper JSON flag file in " + self.raw_flag_dir + " for RUN " + self.run_id))  # + "\nSearching for " + self.run_id + "_" + lane + ".+json"))

        if self.args.splitbarcode_demux:
            # get the barcode names & sequences to add in the JSON flag file
            all_indexes = {}
            for readset in self.readsets[lane]:
                all_indexes[readset.index_name] = readset.index
            with open(raw_flag_file, 'r') as json_fh:
                json_flag_content = json.load(json_fh)
            if self.is_dual_index[lane]:
                json_flag_content['speciesBarcodes'] = dict([(index_name, index_seq['INDEX2']+index_seq['INDEX1']) for index_name, index_seq in all_indexes.items()])
            else:
                json_flag_content['speciesBarcodes'] = dict([(index_name, index_seq['INDEX1']) for index_name, index_seq in all_indexes.items()])
            with open(json_flag_file, 'w') as out_json_fh:
                json.dump(json_flag_content, out_json_fh, indent=4)

        else:
            shutil.copy(raw_flag_file, json_flag_file)

        log.info("JSON FLAG file for lane " + lane + " : " + json_flag_file)
        return json_flag_file

    def generate_lane_sample_sheet(self, lane):
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
        csv_file = os.path.join(self.output_dir, "samplesheet." + lane + ".csv")
        if not os.path.exists(os.path.dirname(csv_file)):
            os.makedirs(os.path.dirname(csv_file))
        writer = csv.DictWriter(
            open(csv_file, 'w'),
            delimiter=str(','),
            fieldnames=csv_headers
        )

        writer.writeheader()
        # barcode validation
        if re.search("B", self.mask[lane]) and not self.args.allow_barcode_collision:
            self.validate_barcodes(lane)
        index_lengths = self.get_smallest_index_length(lane)
        for readset in self.readsets[lane]:
            for idx, readset_index in enumerate(readset.indexes):
                # Barcode sequence should only match with the barcode cycles defined in the mask
                # so we adjust thw lenght of the index sequences accordingly for the "Sample_Barcode" field
                if self.is_dual_index[lane]:
                    sample_barcode = readset_index['INDEX2'][0:index_lengths[0]] + readset_index['INDEX1'][0:index_lengths[1]]
                else:
                    sample_barcode = readset_index['INDEX1'][0:index_lengths[0]]
                if self.last_index < len(sample_barcode):
                    sample_barcode = sample_barcode[0:self.last_index]
                if self.first_index > 1:
                    sample_barcode = sample_barcode[self.first_index-1:]

                readset_index['BARCODE_SEQUENCE'] = sample_barcode
                readset.indexes[idx] = readset_index

                csv_dict = {
                    "Sample_ID": readset_index['SAMPLESHEET_NAME'],
                    "Sample_Name": readset_index['SAMPLESHEET_NAME'] + '_' + readset_index['INDEX_NAME'],
                    "Library_ID": readset_index['LIBRARY'],
                    "Description": readset.name + '_' + readset.library_type + '_' + readset.library_source,
                    "Sample_Barcode": sample_barcode
                }
                writer.writerow(csv_dict)

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

    def generate_basecall_outputs(self, lane):
        basecall_outputs = []
        postprocessing_jobs = []
        unaligned_dir = os.path.join(self.output_dir, "Unaligned." + lane)
        basecall_dir = os.path.join(unaligned_dir, "basecall")

        for readset in self.readsets[lane]:
            readset_r1_outputs = []
            readset_r2_outputs = []
            for index in readset.indexes:
                readset_r1_outputs.extend(
                    os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_" + index['INDEX_NAME'] + "_1.fq.gz"),
                    os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_" + index['INDEX_NAME'] + "_1.fq.fqStat.txt")
                )
                if readset.run_type == "PAIRED_END":
                    readset_r2_outputs.append(
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_" + index['INDEX_NAME'] + "_2.fq.gz"),
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_" + index['INDEX_NAME'] + "_2.fq.fqStat.txt")
                    )

            # If True, then merge the 'Undetermined' reads
            if self.merge_undetermined[lane]:
                readset_r1_outputs.append(
                    os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_undecoded_1.fq.gz"),
                    os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_undecoded_1.fq.fqStat.txt")
                )
                if readset.run_type == "PAIRED_END":
                    readset_r2_outputs.append(
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_undecoded_2.fq.gz"),
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_undecoded_2.fq.fqStat.txt")
                    )
            # Processing R1 fastq outputs :
            #   convert headers from MGI to Illumina format using zcat and awk
            basecall_outputs.extend(readset_r1_outputs)
            postprocessing_jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.dirname(readset.fastq1)),
                        pipe_jobs(
                            [
                                bash.cat(
                                    readset_r1_outputs,
                                    None,
                                    zip=True
                                ),
                                bash.awk(
                                    None,
                                    None,
                                    self.awk_read_1_processing_command()
                                ),
                                bash.gzip(
                                    None,
                                    readset.fastq1
                                )
                            ]
                        )
                    ],
                    name="fastq_convert.R1." + readset.name + "." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
            )

            # Processing "raw R2" fastq (also contains barcode sequences) :
            #   convert headers from MGI to Illumina format
            #   while extracting I1 and I2 to build clean R2, I1 and R2 fastq
            #   using zcat and awk
            if readset.run_type == "PAIRED_END":
                basecall_outputs.extend(readset_r2_outputs)
                outputs = [readset.fastq2, readset.index_fastq1]
                if self.is_dual_index[lane]:
                    outputs.append(readset.index_fastq2)
                postprocessing_jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(os.path.dirname(readset.fastq2)),
                            pipe_jobs(
                                [
                                    bash.cat(
                                        readset_r2_outputs,
                                        None,
                                        zip=True
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        self.awk_read_2_processing_command(readset, lane)
                                    )
                                ]
                            )
                        ],
                        output_dependency=outputs,
                        name="fastq_convert.R2." + readset.name + "." + self.run_id + "." + lane,
                        samples=self.samples[lane]
                    )
                )
        # Process undetermined reads fastq files
        unmatched_R1_fastq = os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_undecoded_1.fq.gz")
        if unmatched_R1_fastq not in basecall_outputs:
            basecall_outputs.append(unmatched_R1_fastq)
        postprocessing_jobs.append(
            pipe_jobs(
                [
                    bash.cat(
                        unmatched_R1_fastq,
                        None,
                        zip=True
                    ),
                    bash.awk(
                        None,
                        None,
                        self.awk_read_1_processing_command()
                    ),
                    bash.gzip(
                        None,
                        os.path.join(unaligned_dir, "Undetermined_S0_L00" + lane + "_R1_001.fastq.gz")
                    )
                ],
                name= "fastq_convert.R1.unmatched." + self.run_id + "." + lane,
                samples=self.samples[lane]
            )
        )
        if self.is_paired_end[lane]:
            unmatched_R2_fastq = os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_undecoded_2.fq.gz")
            if unmatched_R2_fastq not in basecall_outputs:
                basecall_outputs.append(unmatched_R2_fastq)
            outputs = [
                os.path.join(unaligned_dir, "Undetermined_S0_L00" + lane + "_R2_001.fastq.gz"),
                os.path.join(unaligned_dir, "Undetermined_S0_L00" + lane + "_I1_001.fastq.gz")
            ]
            if self.is_dual_index[lane]:
                outputs.append(os.path.join(unaligned_dir, "Undetermined_S0_L00" + lane + "_I2_001.fastq.gz"))
            postprocessing_jobs.append(
                pipe_jobs(
                    [
                        bash.cat(
                            unmatched_R2_fastq,
                            None,
                            zip=True
                        ),
                        bash.awk(
                            None,
                            None,
                            self.awk_read_2_processing_command(None, lane)
                        )
                    ],
                    output_dependency=outputs,
                    name="fastq_convert.R2.unmatched." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
            )

        return basecall_outputs, postprocessing_jobs

    def generate_demuxfastqs_outputs(self, lane):
        demuxfastqs_outputs = []
        postprocessing_jobs = []
        output_dir = os.path.join(self.output_dir, "Unaligned." + lane)

        for readset in self.readsets[lane]:
            readset_r1_outputs = []
            readset_r2_outputs = []

            for index in readset.indexes:
                readset_r1_outputs.append(
                    os.path.join(output_dir, "tmp", index['SAMPLESHEET_NAME']+"-"+index['SAMPLESHEET_NAME']+'_'+index['INDEX_NAME']+"-"+index['BARCODE_SEQUENCE']+"_R1.fastq.gz")
                )
                if readset.run_type == "PAIRED_END":
                    readset_r2_outputs.append(
                        os.path.join(output_dir, "tmp", index['SAMPLESHEET_NAME']+"-"+index['SAMPLESHEET_NAME']+'_'+index['INDEX_NAME']+"-"+index['BARCODE_SEQUENCE']+"_R2.fastq.gz")
                    )

            # If True, then merge the 'Undetermined' reads
            if self.merge_undetermined[lane]:
                readset_r1_outputs.append(
                    os.path.join(output_dir, "tmp", "unmatched_R1.fastq.gz")
                )
                if readset.run_type == "PAIRED_END":
                    readset_r2_outputs.append(
                        os.path.join(output_dir, "tmp", "unmatched_R2.fastq.gz")
                    )

            # Processing R1 fastq outputs :
            #   convert headers from MGI to Illumina format using zcat and awk
            demuxfastqs_outputs.extend(readset_r1_outputs)
            postprocessing_jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.dirname(readset.fastq1)),
                        pipe_jobs(
                            [
                                bash.cat(
                                    readset_r1_outputs,
                                    None,
                                    zip=True
                                ),
                                bash.awk(
                                    None,
                                    None,
                                    self.awk_read_1_processing_command()
                                ),
                                bash.gzip(
                                    None,
                                    readset.fastq1
                                )
                            ]
                        )
                    ],
                    name="fastq_convert.R1." + readset.name + "." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
            )

            # Processing "raw R2" fastq (also contains barcode sequences) :
            #   convert headers from MGI to Illumina format
            #   while extracting I1 and I2 to build clean R2, I1 and R2 fastq
            #   using zcat and awk
            if readset.run_type == "PAIRED_END":
                demuxfastqs_outputs.extend(readset_r2_outputs)
                outputs = [readset.fastq2, readset.index_fastq1]
                if self.is_dual_index[lane]:
                    outputs.append(readset.index_fastq2)
                postprocessing_jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(os.path.dirname(readset.fastq2)),
                            pipe_jobs(
                                [
                                    bash.cat(
                                        readset_r2_outputs,
                                        None,
                                        zip=True
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        self.awk_read_2_processing_command(readset, lane)
                                    )
                                ]
                            )
                        ],
                        output_dependency=outputs,
                        name="fastq_convert.R2." + readset.name + "." + self.run_id + "." + lane,
                        samples=self.samples[lane]
                    )
                )
        # Process undetermined reads fastq files
        unmatched_R1_fastq = os.path.join(output_dir, "tmp", "unmatched_R1.fastq.gz")
        if unmatched_R1_fastq not in demuxfastqs_outputs:
            demuxfastqs_outputs.append(unmatched_R1_fastq)
        postprocessing_jobs.append(
            pipe_jobs(
                [
                    bash.cat(
                        unmatched_R1_fastq,
                        None,
                        zip=True
                    ),
                    bash.awk(
                        None,
                        None,
                        self.awk_read_1_processing_command()
                    ),
                    bash.gzip(
                        None,
                        os.path.join(output_dir, "Undetermined_S0_L00" + lane + "_R1_001.fastq.gz")
                    )
                ],
                name= "fastq_convert.R1.unmatched." + self.run_id + "." + lane,
                samples=self.samples[lane]
            )
        )

        unaligned_i1 = os.path.join(output_dir, "Undetermined_S0_L00" + lane + "_I1_001.fastq.gz")
        unexpected_barcode_counts_i1 = re.sub(".fastq.gz", ".counts.txt", unaligned_i1)
        demuxfastqs_outputs.append(unexpected_barcode_counts_i1)
        postprocessing_jobs.append(
            run_processing_tools.fastq_unexpected_count(
                unaligned_i1,
                unexpected_barcode_counts_i1,
                name="fastq_countbarcodes.I1.unmatched." + self.run_id + "." + lane,
            )
        )

        if self.is_paired_end[lane]:
            unmatched_R2_fastq = os.path.join(output_dir, "tmp", "unmatched_R2.fastq.gz")
            if unmatched_R2_fastq not in demuxfastqs_outputs:
                demuxfastqs_outputs.append(unmatched_R2_fastq)
            outputs = [
                os.path.join(output_dir, "Undetermined_S0_L00" + lane + "_R2_001.fastq.gz"),
                os.path.join(output_dir, "Undetermined_S0_L00" + lane + "_I1_001.fastq.gz")
            ]
            if self.is_dual_index[lane]:
                unaligned_i2 = os.path.join(output_dir, "Undetermined_S0_L00" + lane + "_I2_001.fastq.gz")
                outputs.append(unaligned_i2)
                unexpected_barcode_counts_i2 = re.sub(".fastq.gz", ".counts.txt", unaligned_i2)
                demuxfastqs_outputs.append(unexpected_barcode_counts_i2)
                postprocessing_jobs.append(
                    run_processing_tools.fastq_unexpected_count(
                        unaligned_i2,
                        unexpected_barcode_counts_i2,
                        name="fastq_countbarcodes.I2.unmatched." + self.run_id + "." + lane,
                    )
                )

            postprocessing_jobs.append(
                pipe_jobs(
                    [
                        bash.cat(
                            unmatched_R2_fastq,
                            None,
                            zip=True
                        ),
                        bash.awk(
                            None,
                            None,
                            self.awk_read_2_processing_command(None, lane)
                        )
                    ],
                    output_dependency=outputs,
                    name="fastq_convert.R2.unmatched." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
            )

        return demuxfastqs_outputs, postprocessing_jobs

    def awk_read_1_processing_command(self):
        """
        Returns a string serving as instructions for awk.
        This produces the command to convert the header of R1 fastq file from MGI to Illumina format
        """
        return """-v inst=\"{instrument}\" -v run=\"{run}\" 'match($0, /@({flowcell})L([0-9]{{1}})C([0-9]{{3}})R([0-9]{{3}})([0-9]+):?([ACTGN-]+)?\/([0-9]{{1}})/, head_items) {{
 gsub("^0*", "", head_items[3])
 gsub("^0*", "", head_items[4])
 gsub("^0*", "", head_items[5])
 print "@" inst ":" run ":" head_items[1] ":" head_items[2] ":" head_items[5] ":" head_items[3] ":" head_items[4] " " head_items[7] ":N:0:" head_items[6]
 next
}} 1'""".format(
            instrument=self.instrument,
            run=self.run_number,
            flowcell=self.flowcell_id
        )

    def awk_read_2_processing_command(self, readset, lane):
        """
        Returns a string serving as instructions for awk.
        This produces the command to extract I1 (and I2 if exists) sequence from the R2 fastq,
        creating R2 (without barcode sequences) I1 and I2 fastq files from R2 fastq file.
        This will also convert the headers of the fastqs from MGI Illumina format
        """
        if self.is_dual_index[lane]:
            return """-v inst=\"{instrument}\" -v run=\"{run}\" -v read_len=\"{read_len}\" -v barcode1_len=\"{barcode1_len}\" -v barcode2_len=\"{barcode2_len}\" '{{
 header=$0
 getline seq
 getline sep
 getline qual
 match(header, /@({flowcell})L([0-9]{{1}})C([0-9]{{3}})R([0-9]{{3}})([0-9]+):?([ACTGN-]+)?\/([0-9]{{1}})/, head_items)
 gsub("^0*", "", head_items[3]); gsub("^0*", "", head_items[4]); gsub("^0*", "", head_items[5])
 header="@" inst ":" run ":" head_items[1] ":" head_items[2] ":" head_items[5] ":" head_items[3] ":" head_items[4] " " head_items[7] ":N:0:" head_items[6]
 r2_seq=substr(seq,1,read_len)
 i1_seq=substr(seq,read_len+barcode2_len+1,barcode1_len)
 i2_seq=substr(seq,read_len+1,barcode2_len)
 r2_qual=substr(qual,1,read_len)
 i1_qual=substr(qual,read_len+barcode2_len+1,barcode1_len)
 i2_qual=substr(qual,read_len+1,barcode2_len)
 print header "\\n" r2_seq "\\n" sep "\\n" r2_qual | "gzip > {r2_out}"
 print header "\\n" i1_seq "\\n" sep "\\n" i1_qual | "gzip > {i1_out}"
 print header "\\n" i2_seq "\\n" sep "\\n" i2_qual | "gzip > {i2_out}"
}}'""".format(
                instrument=self.instrument,
                run=self.run_number,
                read_len=self.read2cycles[lane],
                barcode1_len=self.index1cycles[lane],
                barcode2_len=self.index2cycles[lane],
                flowcell=self.flowcell_id,
                r2_out=readset.fastq2 if readset else os.path.join(self.output_dir, "Unaligned." + lane, "Undetermined_S0_L00" + lane + "_R2_001.fastq.gz"),
                i1_out=readset.index_fastq1 if readset else os.path.join(self.output_dir, "Unaligned." + lane, "Undetermined_S0_L00" + lane + "_I1_001.fastq.gz"),
                i2_out=readset.index_fastq2 if readset else os.path.join(self.output_dir, "Unaligned." + lane, "Undetermined_S0_L00" + lane + "_I2_001.fastq.gz")
            )
        else:
            return """-v inst=\"{instrument}\" -v run=\"{run}\" -v read_len=\"{read_len}\" -v barcode_len=\"{barcode_len}\" '{{
 header=$0
 getline seq
 getline sep
 getline qual
 match(header, /@({flowcell})L([0-9]{{1}})C([0-9]{{3}})R([0-9]{{3}})([0-9]+):?([ACTGN]+)?\/([0-9]{{1}})/, head_items)
 gsub("^0*", "", head_items[3]); gsub("^0*", "", head_items[4]); gsub("^0*", "", head_items[5])
 header="@" inst ":" run ":" head_items[1] ":" head_items[2] ":" head_items[5] ":" head_items[3] ":" head_items[4] " " head_items[7] ":N:0:" head_items[6]
 r2_seq=substr(seq,1,read_len)
 i1_seq=substr(seq,read_len+1,barcode_len)
 r2_qual=substr(qual,1,read_len)
 i1_qual=substr(qual,read_len+1,barcode_len)
 print header "\\n" r2_seq "\\n" sep "\\n" r2_qual | "gzip > {r2_out}"
 print header "\\n" i1_seq "\\n" sep "\\n" i1_qual | "gzip > {i1_out}"
}}'""".format(
                instrument=self.instrument,
                run=self.run_number,
                read_len=self.read2cycles[lane],
                barcode_len=self.index1cycles[lane],
                flowcell=self.flowcell_id,
                r2_out=readset.fastq2 if readset else os.path.join(self.output_dir, "Unaligned." + lane, "Undetermined_S0_L00" + lane + "_R2_001.fastq.gz"),
                i1_out=readset.index_fastq1 if readset else os.path.join(self.output_dir, "Unaligned." + lane, "Undetermined_S0_L00" + lane + "_I1_001.fastq.gz")
            )



    def get_smallest_index_length(self, lane):
        """
        Returns a list (for each index read of the run) of the minimum between the number of index cycle on the
        sequencer and all the index lengths.
        """
        run_index_lengths = []
        all_indexes = []
        for readset in self.readsets[lane]:
            all_indexes += readset.index

        if self.is_dual_index[lane]:
            min_sample_index_length = min(len(index['INDEX2']) for index in all_indexes)
            run_index_lengths.append(min(min_sample_index_length, int(self.index2cycles[lane])))

        min_sample_index_length = min(len(index['INDEX1']) for index in all_indexes)
        run_index_lengths.append(min(min_sample_index_length, int(self.index1cycles[lane])))

        return run_index_lengths

    def load_readsets(self, lane):
        """
        Parse the sample sheet and return a list of readsets.
        """

        return parse_mgi_raw_readset_files(
            self.readset_file,
            "PAIRED_END" if self.is_paired_end[lane] else "SINGLE_END",
            self.seqtype,
            self.run_id,
            self.flowcell_id,
            lane,
            int(self.read1cycles[lane]),
            int(self.read2cycles[lane]),
            int(self.index1cycles[lane]),
            int(self.index2cycles[lane]),
            self.output_dir
        )

    def submit_jobs(self):
        super(MGIRunProcessing, self).submit_jobs()

    @property
    def steps(self):
        return [
            [
                self.index,
                self.fastq,
                self.qc_graphs,
                self.fastp,
                self.fastqc,
                self.blast,
                self.align,
                self.picard_mark_duplicates,
                self.metrics,
                self.md5,
                self.report,
                self.copy,
                self.final_notification
            ],
            [
                self.basecall,
                self.fastq,
                self.qc_graphs,
                self.fastp,
                self.fastqc,
                self.blast,
                self.align,
                self.picard_mark_duplicates,
                self.metrics,
                self.md5,
                self.report,
                self.copy,
                self.final_notification
            ]
        ]

def distance(
    str1,
    str2
    ):
    """
    Returns the hamming distance. http://code.activestate.com/recipes/499304-hamming-distance/#c2
    """
    return sum(map(str.__ne__, str1, str2))

if __name__ == '__main__':

    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        MGIRunProcessing(protocol=['g400', 't7'])

