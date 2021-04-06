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
import itertools
import subprocess
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
        self.copy_job_inputs = {}
        self.argparser.add_argument("-r", "--readsets", help="Sample sheet for the MGI run to process (mandatory)", type=file, required=False)
        self.argparser.add_argument("-d", "--run", help="Run directory (mandatory)", required=False, dest="run_dir")
        self.argparser.add_argument("--run-id", help="Run ID. Default is parsed from the run folder", required=False, dest="run_id")
        self.argparser.add_argument("--flowcell-id", help="Flowcell ID. Default is parsed from the run folder", required=False, dest="flowcell_id")
        self.argparser.add_argument("--raw-fastq-prefix", help="Prefix used to search for the raw fastq from the sequencer. Default <FLOWCELL_ID>_<RUN_ID>", required=False, dest="raw_fastq_prefix")
        self.argparser.add_argument("--lane", help="Lane number (to only process the given lane)", type=int, required=False, dest="lane_number")
        self.argparser.add_argument("--demux-fastq", help="Fastq files given by the sequencer are already demultiplexed : NO DEMULTIPLEXING will be performed by the pipeline", action="store_true", required=False, dest="demux_fastq")
        self.argparser.add_argument("-x", help="First index base to use for demultiplexing (inclusive). The index from the sample sheet will be adjusted according to that value.", type=int, required=False, dest="first_index")
        self.argparser.add_argument("-y", help="Last index base to use for demultiplexing (inclusive)", type=int, required=False, dest="last_index")
        self.argparser.add_argument("-m", help="Number of index mistmaches allowed for demultiplexing (default 1). Barcode collisions are always checked.", type=int, required=False, dest="number_of_mismatches")
        self.argparser.add_argument("--allow-barcode-collision", help="Allow barcode collision by not comparing barcode sequences to each other (usually decreases the demultiplexing efficiency).", action="store_true", required=False, dest="allow_barcode_collision")

        super(MGIRunProcessing, self).__init__(protocol)

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = {}
            if not hasattr(self, "_mask"):
                self._mask = {}
            for lane in self.lanes:
                self._readsets[lane] = self.load_readsets(lane)
                self._mask[lane] = self.get_mask(lane)
                self.generate_mgi_lane_sample_sheet(lane)
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
                self._lanes = [lane for lane in list(set([line['Position'].split(":")[0] for line in csv.DictReader(open(self.readset_file, 'rb'), delimiter='\t', quotechar='"')]))]
            for lane in self._lanes:
                self.copy_job_inputs[lane] = []
        return self._lanes

    @property
    def is_paired_end(self):
        if not hasattr(self, "_is_paired_end"):
            self._is_paired_end = {}
            for lane in self.lanes:
                if self.get_read2cycles(lane):
                    self._is_paired_end[lane] = True
                else:
                    self._is_paired_end[lane] = False
        return self._is_paired_end

    @property
    def is_dual_index(self):
        if not hasattr(self, "_is_dual_index"):
            self._is_dual_index = {}
            for lane in self.lanes:
                if self.get_index2cycles(lane) == "0":
                    self._is_dual_index[lane] = False
                else:
                    self._is_dual_index[lane] = True
        return self._is_dual_index

    @property
    def is_demultiplexed(self):
        if not hasattr(self, "_is_demultiplexed"):
            if self.args.demux_fastq:
                self._is_demultiplexed = True
            else:
                self._is_demultiplexed = False
        return self._is_demultiplexed

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
    def is_dual_index(self):
        if not hasattr(self, "_is_dual_index"):
            self._is_dual_index = {}
            for lane in self.lanes:
                if self.get_index2cycles(lane) == "0":
                    self._is_dual_index[lane] = False
                else:
                    self._is_dual_index[lane] = True
        return self._is_dual_index

    @property
    def is_demultiplexed(self):
        if not hasattr(self, "_is_demultiplexed"):
            if self.args.demux_fastq:
                self._is_demultiplexed = True
            else:
                self._is_demultiplexed = False
        return self._is_demultiplexed

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
    def run_counter(self):
        if not hasattr(self, "_run_counter"):
            run_counter = ""
            for lane in self.lanes:
                lane_run_counter = self.get_run_counter(lane)
                if run_counter and run_counter != lane_run_counter:
                    _raise(SanitycheckError("One run (\"" + self.run_id + "\") cannot be defined by more than one run counter !! (\"" + lane_run_counter + "\" vs. \"" +  run_counter + "\")"))
                else:
                    run_counter = lane_run_counter
            self._run_counter = run_counter 
        return self._run_counter

    @property
    def seqtype(self):
        if not hasattr(self, "_seqtype"):
            self._seqtype = self.get_seqtype()
        return self._seqtype

    @property
    def year(self):
        """
        Get year of the from sample sheet
        """
        if not hasattr(self, "_year"):
            dates = set([date for date in list(set([line['Start Date'] for line in csv.DictReader(open(self.readset_file, 'rb'), delimiter='\t', quotechar='"')]))])
            if len(list(dates)) > 1:
                _raise(SanitycheckError("More than one date were found in the sample sheet for the run \"" + self._run_id + "\""))
            else:
                self._year = list(dates)[0].split("-")[0]
        return self._year

    @property
    def date(self):
        """
        Get whole date of the run from sample sheet
        """
        if not hasattr(self, "_date"):
            dates = set([date for date in list(set([line['Start Date'] for line in csv.DictReader(open(self.readset_file, 'rb'), delimiter='\t', quotechar='"')]))])
            if len(list(dates)) > 1:
                _raise(SanitycheckError("More than one date were found in the sample sheet for the run \"" + self._run_id + "\""))
            else:
                date = list(dates)[0].split("-")
                self._date = date[0][-2:] + date[1] + date[2]
        return self._date

    @property
    def run_dir(self):
        if self.args.run_dir:
            return self.args.run_dir
        else:
            _raise(SanitycheckError("Error: missing '-d/--run_dir' option!"))

    @property
    def lane_number(self):
        if self.args.lane_number:
            return self.args.lane_number
        else:
            return None

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
    def lane_number(self):
        if self.args.lane_number:
            return self.args.lane_number
        else:
            return None

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
    def raw_fastq_prefix(self):
        """
        """
        if not hasattr(self, "_raw_fastq_prefix"):
            if self.args.raw_fastq_prefix:
                self._raw_fastq_prefix = self.args.raw_fastq_prefix
            else:
                if self.is_demultiplexed:
                    self._raw_fastq_prefix = self.flowcell_id
                else:
                    self._raw_fastq_prefix = self.flowcell_id + "_" + self.run_id
        return self._raw_fastq_prefix

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
        if self.args.readsets:
            return self.args.readsets.name
        else:
            _raise(SanitycheckError("Error: missing '-r/--readsets' argument !"))

    @property
    def bioinfo_files(self):
        if not hasattr(self, "_bioinfo_files"):
            self._bioinfo_files = {}
            for lane in self.lanes:
                if os.path.exists(os.path.join(self.output_dir, "L0" + lane, "Unaligned." + lane, "raw_fastq", "BioInfo.csv")):
                    self._bioinfo_files[lane] = os.path.join(self.output_dir, "L0" + lane, "Unaligned." + lane, "raw_fastq", "BioInfo.csv")
                else:
                    rundir = self.run_dir
                    if "(no_barcode)" in rundir:
                        rundir = rundir.split("(")[0]
                    bioinfo_file = os.path.join(rundir, "L0" + lane, "BioInfo.csv")
                    if lane in self._bioinfo_files and self._bioinfo_files[lane] != bioinfo_file:
                        _raise(SanitycheckError("More than one Bioinfo.csv found for lane '" + lane + "' : " + self._bioinfo_files[lane] + ", " + bioinfo_file))
                    else:
                        self._bioinfo_files[lane] = bioinfo_file
        return self._bioinfo_files

    @property
    def report_hash(self):
        if not hasattr(self, "_report_hash"):
            self._report_hash = {}
            for lane in self.lanes:
                self._report_hash[lane] = {
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
                    'index': {},
                    'qc' : {},
                    'blast' : {},
                    'align' : {},
                    'mark_dup' : {}
                }
        return self._report_inputs

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
            lane_jobs = []

            # First we copy all the lane folder into the raw_fastq folder
            unaligned_dir = os.path.join(self.output_dir, "L0" + lane, "Unaligned." + lane)
            raw_fastq_dir = os.path.join(unaligned_dir, "raw_fastq")
            copy_done_file = os.path.join(raw_fastq_dir, "copy_done.Success")

            raw_name_prefix = self.raw_fastq_prefix +  "_L0" + lane
            copy_job_output_dependency = [
                os.path.join(raw_fastq_dir, raw_name_prefix + ".summaryReport.html"),
                os.path.join(raw_fastq_dir, raw_name_prefix + ".heatmapReport.html"),
                os.path.join(raw_fastq_dir, "summaryTable.csv")
            ]
            rename_job = None
    
            # If demultiplexing was done on the sequencer i.e. MGI adpaters used...
            if self.is_demultiplexed:
                rename_jobs = []
                copy_job_output_dependency.extend([
                    os.path.join(raw_fastq_dir, "BarcodeStat.txt"),
                    os.path.join(raw_fastq_dir, "SequenceStat.txt")
                ])
    
                # ...then do the renamings from the raw_fastq folder to the Unaligned folder
                for readset in self.readsets[lane]:
                    output_dir = os.path.join(self.output_dir, "L0" + lane, "Unaligned." + lane, 'Project_' + readset.project_id, 'Sample_' + readset.name)
                    m = re.search("\D+(?P<index>\d+)", readset.index_name)
                    if m:
                        index_number = m.group('index').lstrip("0")
                    else:
                        _raise(SanitycheckError("MGI Index error : unable to parse barcode # in : " + readset.index_name))
    
                    copy_job_output_dependency.extend([
                        os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + lane + "_" + index_number + "_1.fq.gz"),
                        os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + lane + "_" + index_number + "_1.fq.fqStat.txt"),
                        os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + lane + "_" + index_number + ".report.html")
                    ])
                    rename_jobs.append(
                        concat_jobs([
                            bash.mkdir(output_dir),
                            bash.cp(
                                os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + lane + "_" + index_number + "_1.fq.gz"),
                                readset.fastq1
                            ),
                            bash.cp(
                                os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + lane + "_" + index_number + "_1.fq.fqStat.txt"),
                                re.sub("gz", "fqStat.txt", readset.fastq1)
                            ),
                            bash.cp(
                                os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + lane + "_" + index_number + ".report.html"),
                                re.sub("_R1_001.fastq.gz", ".report.html", readset.fastq1)
                            )
                        ])
                    )
    
                    if readset.run_type == "PAIRED_END":
                        copy_job_output_dependency.extend([
                            os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + lane + "_" + index_number + "_2.fq.gz"),
                            os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + lane + "_" + index_number + "_2.fq.fqStat.txt")
                        ])
                        rename_jobs.append(
                            concat_jobs([
                                bash.cp(
                                    os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + lane + "_" + index_number + "_2.fq.gz"),
                                    readset.fastq2,
                                ),
                                bash.cp(
                                    os.path.join(raw_fastq_dir, readset.flow_cell + "_L0" + lane + "_" + index_number + "_2.fq.fqStat.txt"),
                                    re.sub("gz", "fqStat.txt", readset.fastq2),
                                )
                            ])
                        )
    
                rename_job = concat_jobs(rename_jobs)
                rename_job.output_files.append(copy_done_file)
                rename_job.name = "index.rename." + self.run_id + "." + lane
                rename_job.samples = self.samples[lane]
    
            elif "(no_barcode)" in self.run_dir:
                if len(self.readsets[lane]) > 1:
                    err_msg = "LANE SETTING ERROR :\n"
                    err_msg += "Unable to demultiplex " + str(len(self.readsets[lane])) + " samples : No barcode in fastq files...\n(in "
                    err_msg += self.run_dir + ")"
                    _raise(SanitycheckError(err_msg))
                readset = self.readsets[lane][0]
                rename_job = concat_jobs(
                    [
                        bash.mkdir(os.path.dirname(readset.fastq1)),
                        bash.cp(
                            os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + lane + "_read_1.fq.gz"),
                            readset.fastq1,
                        ),
                        bash.cp(
                            os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + lane + "_read_1.fq.fqStat.txt"),
                            re.sub("gz", "fqStat.txt", readset.fastq1),
                        ),
                        bash.cp(
                            os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + lane + "_read.report.html"),
                            re.sub("_R1_001.fastq.gz", "_read.report.html", readset.fastq1),
                        )
                    ]
                )
                if readset.run_type == "PAIRED_END":
                    rename_job = concat_jobs(
                        [
                            rename_job,
                            bash.cp(
                                os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + lane + "_read_2.fq.gz"),
                                readset.fastq2,
                            ),
                            bash.cp(
                                os.path.join(raw_fastq_dir, self.flowcell_id + "_L0" + lane + "_read_2.fq.fqStat.txt"),
                                re.sub("gz", "fqStat.txt", readset.fastq2),
                            )
                        ]
                    )
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
                # Handle corner cases for BioInfo.csv
                if "(no_barcode)" in self.run_dir:
                    copy_job = concat_jobs(
                        [
                            copy_job,
                            bash.cp(
                                self.bioinfo_files[lane],
                                raw_fastq_dir
                            )
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
    
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)
        return jobs

    def fastq(self):
        """
        *** In the future, may generate the fastq files from the raw CAL files. ***
        Perform demultplexing of the reads with fgbio DemuxFastqs
        """
        jobs = []

        if self.is_demultiplexed:
            log.info("Demultiplexing done on the sequencer... Skipping fastq step...")

        elif "(no_barcode)" in self.run_dir:
            log.info("No barcode in fastq files... Demultiplexing is undoable... Skipping fastq step...")

        else:
            log.info("No demultiplexing done on the sequencer... Processing fastq step...")

            for lane in self.lanes:
                lane_jobs = []
                jobs_to_throttle = []

                input_fastq_dir = os.path.join(self.output_dir, "L0" + lane, "Unaligned." + lane, "raw_fastq")
                raw_name_prefix = self.raw_fastq_prefix +  "_L0" + lane
                input1 = os.path.join(input_fastq_dir, raw_name_prefix + "_read_1.fq.gz")
                input2 = os.path.join(input_fastq_dir, raw_name_prefix + "_read_2.fq.gz")

                demuxfastqs_outputs, postprocessing_jobs = self.generate_demuxfastqs_outputs(lane)

                tmp_output_dir = os.path.dirname(demuxfastqs_outputs[0])    
                tmp_metrics_file = os.path.join(tmp_output_dir, self.run_id + "." + lane + ".DemuxFastqs.metrics.txt")
                metrics_file = os.path.join(self.output_dir, "L0" + lane,  "Unaligned." + lane, self.run_id + "." + lane + ".DemuxFastqs.metrics.txt")
                demuxfastqs_outputs.append(metrics_file)

                demultiplex_job = run_processing_tools.demux_fastqs(
                    os.path.join(self.output_dir, "L0" + lane, "samplesheet." + lane + ".csv"),
                    self.number_of_mismatches,
                    self.mask[lane],
                    demuxfastqs_outputs,
                    tmp_metrics_file,
                    input1,
                    input2
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
                        name="fastq.demultiplex." + self.run_id + "." + lane,
                        samples=self.samples[lane],
                        input_dependency=demultiplex_job.input_files,
                        output_dependency=demuxfastqs_outputs
                    )
                )
    
                if postprocessing_jobs:
                    jobs_to_throttle.extend(postprocessing_jobs)
    
                for readset in self.readsets[lane]:
                    self.report_inputs[lane]['index'][readset.name] = [metrics_file]
   
                lane_jobs.extend(self.throttle_jobs(jobs_to_throttle)) 
                self.add_copy_job_inputs(lane_jobs, lane)
                jobs.extend(lane_jobs)
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
                    name="qc." + readset.name + ".qc." + self.run_id + "." + lane,
                    samples=[readset.sample]
                ))
                self.report_inputs[lane]['qc'][readset.name] = os.path.join(output_dir, "mpsQC_" + region_name + "_stats.xml")
   
                if not self.is_demultiplexed:
                    # Also process the fastq of the indexes
                        lane_jobs.append(
                            concat_jobs([
                            bash.mkdir(output_dir + "_index"),
                            bvatools.readsqc(
                                readset.index_fastq1,
                                readset.index_fastq2 if self.is_dual_index[lane] else None,
                                "FASTQ",
                                region_name,
                                output_dir + "_index"
                        )],
                        name="qc." + readset.name + ".qc_index." + self.run_id + "." + lane,
                        samples=[readset.sample]
                    ))
    
            self.add_to_report_hash("qc_graphs", lane, lane_jobs)
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
                        os.path.join(os.path.dirname(readset.fastq1), "fastqc.R1", re.sub(".fastq.gz", "_fastqc", os.path.basename(readset.fastq1)), re.sub(".fastq.gz", "_fastqc", "fastqc_data.txt")) 
                    ]
                }
                if readset.run_type == "PAIRED_END":
                    input_dict[readset.fastq2] = [
                        os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.fastq2))),
                        os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.fastq2))),
                        os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc", os.path.basename(readset.fastq2)), re.sub(".fastq.gz", "_fastqc", "fastqc_data.txt"))
                    ]
                if not self.is_demultiplexed:
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
                        if self.is_demultiplexed:
                            self.report_inputs[lane]['index'][readset.name] = [outputs[0]]
                        else:
                            self.report_inputs[lane]['index'][readset.name].append(outputs[0])
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

        for lane in self.lanes:
            lane_jobs = []

            if is_nb_blast_per_lane:
                nb_blast_to_do = int(nb_blast_to_do) // len(self.readsets[lane])

            nb_blast_to_do = max(1, nb_blast_to_do)

            for readset in self.readsets[lane]:
                output_prefix = os.path.join(
                    self.output_dir,
                    "L0" + lane,
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
                    name="blast." + readset.name + ".blast." + self.run_id + "." + lane,
                    samples = [readset.sample]
                )
                lane_jobs.append(job)
    
            self.add_to_report_hash("blast", lane, lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)
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
    
            self.add_to_report_hash("picard_mark_duplicates", lane, lane_jobs)
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
    
            self.add_to_report_hash("metrics", lane, lane_jobs)
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
                lane_job = concat_jobs(
                    lane_jobs,
                    name="md5." + self.run_id + "." + lane,
                    samples=self.readsets[lane]
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

            report_dir = os.path.join(self.output_dir, "L0" + lane, "report")
            sample_report_dir  = os.path.join(report_dir, "sample_json")
    
            run_validation_inputs = []
    
            # Add barcodes info to the report_hash
            self.report_hash[lane]["barcodes"] = dict([(readset.name, readset.indexes) for readset in self.readsets[lane]])
    
            general_information_file = os.path.join(self.output_dir, "L0" + lane, self.run_id + "." + lane + ".general_information.json")
            with open(general_information_file, 'w') as out_json:
                json.dump(self.report_hash[lane], out_json, indent=4)
    
            # Build JSON report for each sample
            for readset in self.readsets[lane]:
                # Build JSON report for each sample
                output_file = os.path.join(sample_report_dir, readset.name + ".report.json")
                lane_jobs.append(
                    concat_jobs([
                        bash.mkdir(sample_report_dir),
                        tools.run_validation_sample_report(
                            readset,
                            self.report_inputs[lane],
                            output_file,
                            general_information_file
                        )],
                        name="sample_report." + readset.name + "." + self.run_id + "." + lane
                    )
                )
                run_validation_inputs.append(output_file)
    
            run_validation_report_json = os.path.join(report_dir, self.run_id + "." + lane + ".run_validation_report.json")
            # Aggregate sample reports to build the run validation report
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
    
            # Copy fastqc HTML files into the report folder
            copy_fastqc_jobs = [
                bash.mkdir(os.path.join(self.output_dir, "L0" + lane, "report", "fastqc"))
            ]
            for readset in self.readsets[lane]:
                copy_fastqc_jobs.append(
                    concat_jobs([
                        bash.cp(
                            os.path.join(os.path.dirname(readset.fastq1), "fastqc.R1", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.fastq1))),
                            os.path.join(self.output_dir, "L0" + lane, "report", "fastqc/")
                        ),
                        bash.cp(
                            os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.fastq2))),
                            os.path.join(self.output_dir, "L0" + lane, "report", "fastqc/")
                        ) if readset.run_type == "PAIRED_END" else None
                    ]
                )
            )
            copy_fastqc_jobs[0].output_files.append(os.path.join(self.output_dir, "L0" + lane, "report", "fastqc"))
    
            # Copy the MGI summary HTML report into the report folder
            raw_name_prefix = self.raw_fastq_prefix +  "_L0" + lane
            summary_report_html = os.path.join(self.output_dir, "L0" + lane, "Unaligned." + lane, "raw_fastq", raw_name_prefix + ".summaryReport.html")
            lane_jobs.append(
                concat_jobs(
                    copy_fastqc_jobs + [
                    bash.cp(
                        summary_report_html,
                        os.path.join(self.output_dir, "L0" + lane, "report")
                    )],
                    name="report.copy." + self.run_id + "." + lane
                )
            )
    
            # Create a report zip archive containing all the above...
            zip_output = os.path.join(self.output_dir, "L0" + lane, "report", self.run_id + "_" + self.flowcell_id + "_L00" + lane + "_report.zip")
            zip_job = bash.zip(
                [
                    run_validation_report_json,
                    os.path.join(self.output_dir, "L0" + lane, "report", "fastqc"),
                    summary_report_html
                ],
                zip_output,
                recursive=True
            )
            zip_job.name = "report.zip." + self.run_id + "." + lane
            zip_job.samples = self.samples[lane]
            lane_jobs.append(zip_job)

            self.add_copy_job_inputs(lane_jobs, lane)

            jobs.extend(lane_jobs)

        return self.throttle_jobs(jobs)

    def copy(self):
        """
        Copy the whole processing foler to where they can be serve or loaded into a LIMS
        """
        jobs_to_concat = []

        full_destination_folder = os.path.join(
            config.param("copy", "destination_folder", type="dirpath"),
            self.seqtype,
            self.year,
            self.date + "_" + self.instrument + "_" + self.run_counter + "_" + self.flowcell_position + self.flowcell_id + "_" + self.sequencer_run_id + "-" + self.seqtype
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
                for readset in [readset for readset in self.readsets if readset.bam]:
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
                source=os.path.join(self.output_dir, "L0" +lane),
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

            # transfer reports to data-hub
            transfer_command_report_zip = config.param('copy', 'transfer_report_command', required=False).format(
                year=self.year,
                file=os.path.join(self.output_dir, "L0" + lane, "report", self.run_id + "_" + self.flowcell_id + "_L00" + lane + "_report.zip")
            )

            jobs_to_concat.append(concat_jobs(
                [
                    Job(
                        command=transfer_command_report_zip
                    ),
                    bash.touch(report_transfer_output)
                ],
                input_dependency=[os.path.join(self.output_dir, "L0" + lane, "report", self.run_id + "_" + self.flowcell_id + "_L00" + lane + "_report.zip")],
                output_dependency=[report_transfer_output]
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
            full_destination_folder = os.path.join(
                config.param("copy", "destination_folder", type="dirpath"),
                self.seqtype,
                self.year,
                self.date + "_" + self.instrument + "_" + self.run_counter + "_" + self.flowcell_position + self.flowcell_id + "_" + self.sequencer_run_id + "-" + self.seqtype
            )
            inputs.append(os.path.join(
                full_destination_folder,
                "copyCompleted." + lane + ".out"
            ))

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
                "input_files" : [os.path.relpath(input_file, os.path.join(self.output_dir, "L0" + lane, "report")) for input_file in job.input_files],
                "output_files" : [os.path.relpath(output_file, os.path.join(self.output_dir, "L0" + lane, "report")) for output_file in job.output_files],
                "command" : job.command,
                "modules" : job.modules
            } for job in jobs]
        }
        self.report_hash[lane]["steps"].append(report_hash)

    def get_read1cycles(self, lane):
        """
        Parse the BioInfo.csv file for the number of read 1 cycles in the run
        """
        bioinfo_csv = csv.reader(open(self.bioinfo_files[lane], 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Read1 Cycles":
                return row[1]
        _raise(SanitycheckError("Could not get Read 1 Cycles from " + self.bioinfo_files[lane]))

    def get_read2cycles(self, lane):
        """
        Parse the BioInfo.csv file for the number of read 2 cycles in the run
        """
        bioinfo_csv = csv.reader(open(self.bioinfo_files[lane], 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Read2 Cycles":
                return row[1]
        _raise(SanitycheckError("Could not get Read 2 Cycles from " + self.bioinfo_files[lane]))

    def get_index1cycles(self, lane):
        """
        Parse the BioInfo.csv file for the number of index 1 cycles in the run
        """
        bioinfo_csv = csv.reader(open(self.bioinfo_files[lane], 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Barcode":
                return row[1]
        _raise(SanitycheckError("Could not get Index 1 Cycles from " + self.bioinfo_files[lane]))

    def get_index2cycles(self, lane):
        """
        Parse the BioInfo.csv file for the number of index 2 cycles in the run
        """
        bioinfo_csv = csv.reader(open(self.bioinfo_files[lane], 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Dual Barcode":
                if row[1] == '':
                    return "0"
                else:
                    return row[1]
        _raise(SanitycheckError("Could not get Index 2 Cycles from " + self.bioinfo_files[lane]))

    def get_instrument(self, lane):
        """
        Parse the BioInfo.csv file for the instrument name the run has been running on
        """
        bioinfo_csv = csv.reader(open(self.bioinfo_files[lane], 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Machine ID":
                return row[1]
        _raise(SanitycheckError("Could not find intrument from " + self.bioinfo_files[lane]))

    def get_flowcell_position(self, lane):
        """
        Parse the BioInfo.csv file for flowcell position of the run
        """
        bioinfo_csv = csv.reader(open(self.bioinfo_files[lane], 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Flowcell Pos":
                return row[1]
        _raise(SanitycheckError("Could not find flowcell position from " + self.bioinfo_files[lane]))

    def get_sequencer_run_id(self, lane):
        """
        Parse the BioInfo.csv file for the ID of the run given by the sequencer
        """
        bioinfo_csv = csv.reader(open(self.bioinfo_files[lane], 'rb'))
        for row in bioinfo_csv:
            if row[0] == "DNB ID":
                # dnb_id format looks like : 10074MG01B_Lane4
                # where run_id is : 10074MG01B
                return row[1].split("_")[0]
        _raise(SanitycheckError("Could not find DNB ID from " + self.bioinfo_files[lane]))

    def get_run_counter(self, lane):
        """
        Parse the BioInfo.csv file for the ID of the run to extract the run counter
        """
        bioinfo_csv = csv.reader(open(self.bioinfo_files[lane], 'rb'))
        for row in bioinfo_csv:
            if row[0] == "DNB ID":
                # dnb_id format looks like : 10074MG01B_Lane4
                # where run_id is : 10074MG01B
                # and run counter is : 10074
                run_counter = row[1].split("_")[0][:-5]
                # sometimes, format is misleading : 1074
                # so we correct it to : 10074
                while len(run_counter) < 5:
                    run_counter[:1] + "0" +  run_counter[1:]
                return run_counter
        _raise(SanitycheckError("Could not find intrument from " + self.bioinfo_files[lane]))

    def get_seqtype(self):
        """
        Determine which kind of sequencing (iseq, miseq, novaseq, hiseqx, hiseq4000, hiseq2500 or dnbseqg400) was performed,
        depending on the instrument used for the run
        """

        instrument = self.instrument
        instrument_file = config.param('DEFAULT', 'instrument_list_file', type='filepath', required=False)
        if not (instrument_file and os.path.isfile(instrument_file)):
            instrument_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'instrument_list.csv')

        return subprocess.check_output("grep -m1 '"+instrument+"' %s | awk -F',' '{print $3}'" % instrument_file, shell=True).strip()

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
    
        index_cycles = [int(self.get_index1cycles(lane))]
        if self.is_dual_index[lane]:
            index_cycles.insert(0, int(self.get_index2cycles(lane)))
    
        mask = self.get_read1cycles(lane) + 'T ' + self.get_read2cycles(lane) + 'T'
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

    def generate_mgi_lane_sample_sheet(self, lane):
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
            "L0" + lane,
            "samplesheet." + lane + ".csv"
        )
        if not os.path.exists(os.path.dirname(csv_file)):
            os.makedirs(os.path.dirname(csv_file))
        writer = csv.DictWriter(
            open(csv_file, 'wb'),
            delimiter=str(','),
            fieldnames=csv_headers
        )

        writer.writeheader()
        # barcode validation
        if re.search("B", self.mask[lane]) and not self.args.allow_barcode_collision:
            self.validate_barcodes(lane)
        index_lengths = self.get_smallest_index_length(lane)
        for readset in self.readsets[lane]:
            for readset_index in readset.indexes:
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

                csv_dict = {
                    "Sample_ID": readset_index['SAMPLESHEET_NAME'],
                    "Sample_Name": readset_index['SAMPLESHEET_NAME'] + '_' + readset_index['INDEX_NAME'],
                    "Library_ID": readset_index['LIBRARY'],
                    "Description": readset.name + '_' + readset.library_type + '_' + readset.library_source,
                    "Sample_Barcode": sample_barcode
                }
                writer.writerow(csv_dict)

    def generate_demuxfastqs_outputs(self, lane):
        demuxfastqs_outputs = []
        postprocessing_jobs = []
        output_dir = os.path.join(self.output_dir, "L0" + lane, "Unaligned." + lane)

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
                    name="fastq.convert_R1." + readset.name + "." + self.run_id + "." + lane,
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
                    outputs.extend(readset.index_fastq2)
                postprocessing_jobs.append(
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
                        ],
                        output_dependency=outputs,
                        name="fastq.convert_R2." + readset.name + "." + self.run_id + "." + lane,
                        samples=self.samples[lane]
                    )
                )

        return demuxfastqs_outputs, postprocessing_jobs

    def awk_read_1_processing_command(self):
        """
        Returns a string serving as instructions for awk.
        This produces the command to convert the header of R1 fastq file from MGI to Illumina format
        """
        return """-v inst=\"{instrument}\" -v run=\"{run}\" 'match($0, /@(V[0-9]+)L([0-9]{{1}})C([0-9]{{3}})R([0-9]{{3}})([0-9]+):([ACTGN]+)\/([0-9]{{1}})/, head_items) {{
 gsub("^0*", "", head_items[3])
 gsub("^0*", "", head_items[4])
 gsub("^0*", "", head_items[5])
 print "@" inst ":" run ":" head_items[1] ":" head_items[2] ":" head_items[5] ":" head_items[3] ":" head_items[4] " " head_items[7] ":N:0:" head_items[6]
 next
}} 1'""".format(
            instrument=self.instrument,
            run=self.run_counter
        )

    def awk_read_2_processing_command(self, readset, lane):
        """
        Returns a string serving as instructions for awk.
        This produces the command to extract I1 (and I2 if exists) sequence from the R2 fastq,
        creating R2 (without barcode sequences) I1 and I2 fastq files from R2 fastq file.
        This will also convert the headers of the fastqs from MGI Illumina format
        """
        if self.is_dual_index[lane]:
            return """-v inst=\"{{instrument}}\" -v run=\"{{run}}\" -v read_len=\"{{read_len}}\" -v barcode_len=\"{{barcode_len}}\" '{{
 header=$0
 getline seq
 getline sep
 getline qual
 match(header, /@(V[0-9]+)L([0-9]{{1}})C([0-9]{{3}})R([0-9]{{3}})([0-9]+):([ACTGN]+)\/([0-9]{{1}})/, head_items)
 gsub("^0*", "", head_items[3]); gsub("^0*", "", head_items[4]); gsub("^0*", "", head_items[5])
 header="@" inst ":" run ":" head_items[1] ":" head_items[2] ":" head_items[5] ":" head_items[3] ":" head_items[4] " " head_items[7] ":N:0:" head_items[6]
 b1_head="@" inst ":" run ":" head_items[1] ":" head_items[2] ":" head_items[5] ":" head_items[3] ":" head_items[4] " " head_items[7] ":N:0:1"
 r2_seq=substr(seq,1,read_len)
 i1_seq=substr(seq,read_len+barcode_len+1,barcode_len)
 i2_seq=substr(seq,read_len+1,barcode_len)
 r2_qual=substr(qual,1,read_len)
 i1_qual=substr(qual,read_len+barcode_len+1,barcode_len)
 i2_qual=substr(qual,read_len+1,barcode_len)
 print header "\n" r2_seq "\n" sep "\n" r2_qual | "gzip > {{r2_out}{
 print b1_head "\n" i1_seq "\n" sep "\n" i1_qual | "gzip > {{i1_out}}
 print header "\n" i2_seq "\n" sep "\n" i2_qual | "gzip > {{i2_out}}
}}'
""".format(
                intrument=self.instrument,
                run=self.run_counter,
                read_len=self.get_read2cycles(lane),
                barcode_len=self.get_index2cycles(lane),
                r2_out=readset.fastq2,
                i1_out=readset.index_fastq1,
                i2_out=readset.index_fastq2
            )
        else:
            return """-v inst=\"{{instrument}}\" -v run=\"{{run}}\" -v read_len=\"{{read_len}}\" -v barcode_len=\"{{barcode_len}}\" '{{
 header=$0
 getline seq
 getline sep
 getline qual
 match(header, /@(V[0-9]+)L([0-9]{{1}})C([0-9]{{3}})R([0-9]{{3}})([0-9]+):([ACTGN]+)\/([0-9]{{1}})/, head_items)
 gsub("^0*", "", head_items[3]); gsub("^0*", "", head_items[4]); gsub("^0*", "", head_items[5])
 header="@" inst ":" run ":" head_items[1] ":" head_items[2] ":" head_items[5] ":" head_items[3] ":" head_items[4] " " head_items[7] ":N:0:" head_items[6]
 b1_head="@" inst ":" run ":" head_items[1] ":" head_items[2] ":" head_items[5] ":" head_items[3] ":" head_items[4] " " head_items[7] ":N:0:1"
 r2_seq=substr(seq,1,read_len)
 i1_seq=substr(seq,read_len+1,barcode_len)
 r2_qual=substr(qual,1,read_len)
 i1_qual=substr(qual,read_len+1,barcode_len)
 print header "\n" r2_seq "\n" sep "\n" r2_qual | "gzip > {{r2_out}}
 print b1_head "\n" i1_seq "\n" sep "\n" i1_qual | "gzip > {{i1_out}}
}}'
""".format(
                intrument=self.instrument,
                run=self.run_counter,
                read_len=self.get_read2cycles(lane),
                barcode_len=self.get_index2cycles(lane),
                r2_out=readset.fastq2,
                i1_out=readset.index_fastq1
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
            run_index_lengths.append(min(min_sample_index_length, int(self.get_index2cycles(lane))))

        min_sample_index_length = min(len(index['INDEX1']) for index in all_indexes)
        run_index_lengths.append(min(min_sample_index_length, int(self.get_index1cycles(lane))))

        return run_index_lengths

    def load_readsets(self, lane):
        """
        Parse the sample sheet and return a list of readsets.
        """

        return parse_mgi_raw_readset_files(
            self.readset_file,
            self.bioinfo_files[lane],
            "PAIRED_END" if self.is_paired_end[lane] else "SINGLE_END",
            self.seqtype,
            self.flowcell_id,
            lane,
            os.path.join(self.output_dir, "L0" + lane)
        )

    def submit_jobs(self):
        super(MGIRunProcessing, self).submit_jobs()

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

    @property
    def steps(self):
        return [
            self.index,
#            self.fastq,
            self.demuxfastq,
            self.fastq,
            self.qc_graphs,
            self.fastqc,
            self.blast,
            self.align,
            self.picard_mark_duplicates,
            self.metrics,
#            self.md5,
            self.report,
            self.copy,
            self.final_notification
        ]

def distance(
    str1,
    str2
    ):
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

