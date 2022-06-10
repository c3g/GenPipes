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
import csv
import errno
import json
import logging
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as Xml
from Bio.Seq import Seq
from collections import Counter, OrderedDict

# Append genpipes directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs

from bfx.readset import parse_illumina_raw_readset_files, parse_mgi_raw_readset_files
from bfx import bvatools
from bfx import picard
from bfx import fastp
from bfx import fastqc
from bfx import tools
from bfx import run_processing_tools
from bfx import bash_cmd as bash
import utils

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

class RunProcessing(common.MUGQICPipeline):
    """
    MGI Run Processing Pipeline
    ================================

    The standard Run Processing pipeline handles both Illumina and MGI sequencing technologies.
    It uses the Illumina bcl2fastq software to convert and demultiplex Illumina base call files to fastq files.
    In the case of MGI run processing, it uses fastq files produced by the MGI-G400 sequencer, or MGI-T7 base call files,
    then does demultiplexing. Finally, the pipeline runs some QCs on the raw data, on the fastq and on the alignment.

    Sample Sheets
    -------------

    The pipeline uses one input sample sheet, a tsv file having the following columns:

    - ProjectLUID
    - ProjectName
    - ContainerName
    - Position
    - Index
    - LibraryLUID
    - LibraryProcess
    - SampleLUID
    - SampleName
    - Reference
    - Start Date
    - Sample Tag
    - Target Cells
    - Species
    - UDF/Genome Size (Mb)
    - Gender
    - Pool Fraction
    - Capture Type
    - Capture Name
    - Capture REF_BED
    - Library Size
    - Library Kit Name
    - Capture Kit Type
    - Capture Bait Version
    - ChIP-Seq Mark

    Example:
    ProjectLUID	ProjectName	ContainerName	Position	Index	LibraryLUID	LibraryProcess	SampleName	Reference	Start Date	Sample Tag	Target Cells	Species	UDF/Genome Size (Mb)	Gender	Pool Fraction	Capture Type	Capture Name	Capture REF_BED	Library Size	Library Kit Name	Capture Kit Type	Capture Bait Version	ChIP-Seq Mark
    AUL208	NA_Control	V300096783	4:1	MGI09_A02_Barcode_41	2-1981727	RNASeq MGI	RNA_GM12878_1	N/A	2021-12-09	N/A	N/A	Eukaryota:Homo sapiens (Taxon ID:9606)	3257.32	F	0.125	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
    AUL208	NA_Control	V300096783	4:1	MGI13_E02_Barcode_45	2-1981728	RNASeq MGI	RNA_GM12878_1	N/A	2021-12-09	N/A	N/A	Eukaryota:Homo sapiens (Taxon ID:9606)	3257.32	F	0.125	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
    AUL208	NA_Control	V300096783	2:1	MGI10_B02_Barcode_42	2-1981725	RNASeq MGI	RNA_Mother_HG004_GM24143_1	N/A	2021-12-09	N/A	N/A	Eukaryota:Homo sapiens (Taxon ID:9606)	3257.32	F	0.125	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
    AUL208	NA_Control	V300096783	2:1	MGI14_F02_Barcode_46	2-1981726	RNASeq MGI	RNA_Mother_HG004_GM24143_1	N/A	2021-12-09	N/A	N/A	Eukaryota:Homo sapiens (Taxon ID:9606)	3257.32	F	0.125	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
    AUL208	NA_Control	V300096783	4:1	MGI10_B02_Barcode_42	2-1981725	RNASeq MGI	RNA_Mother_HG004_GM24143_1	N/A	2021-12-09	N/A	N/A	Eukaryota:Homo sapiens (Taxon ID:9606)	3257.32	F	0.125	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
    """

    def __init__(self, protocol=None):
        self._protocol=protocol
        self.copy_job_inputs = {}

        self.argparser.add_argument("-t", "--type", help = "Sequencing technology : Illumina, MGI G400 or MGI T7 (mandatory)", choices=['illumina', 'mgig400', 'mgit7'], required=False)
        self.argparser.add_argument("-r", "--readsets", help="Sample sheet for the MGI run to process (mandatory)", type=argparse.FileType('r'), required=False)
        self.argparser.add_argument("-d", "--run", help="Run directory (mandatory)", required=False, dest="run_dir")
        self.argparser.add_argument("--run-id", help="Run ID. Default is parsed from the run folder", required=False, dest="run_id")
        self.argparser.add_argument("-f", "--flag", help="T7 flag files directory (mandatory for MGI T7 runs)", type=pathlib.Path, dest="raw_flag_dir", required=False)
        self.argparser.add_argument("--splitbarcode-demux", help="demultiplexing done while basecalling with MGI splitBarcode  (only affect MGI G400 or T7 runs)", action="store_true", required=False, dest="splitbarcode_demux")
        self.argparser.add_argument("--lane", help="Lane number (to only process the given lane)", type=int, required=False, dest="lane_number")
        self.argparser.add_argument("-x", help="First index base to use for demultiplexing (inclusive). The index from the sample sheet will be adjusted according to that value.", type=int, required=False, dest="first_index")
        self.argparser.add_argument("-y", help="Last index base to use for demultiplexing (inclusive)", type=int, required=False, dest="last_index")
        self.argparser.add_argument("-m", help="Number of index mistmaches allowed for demultiplexing (default 1). Barcode collisions are always checked.", type=int, required=False, dest="number_of_mismatches")
        self.argparser.add_argument("--allow-barcode-collision", help="Allow barcode collision by not comparing barcode sequences to each other (usually decreases the demultiplexing efficiency).", action="store_true", required=False, dest="allow_barcode_collision")

        args = sys.argv[1:]
        typearg = ""
        flag = False
        splitbarcode = False
        for i, arg in enumerate(args):
            if arg in ['-t', '--type']:
                typearg = args[i+1]
            if arg in ['-f', '--flag']:
                flag = True
            if arg == '--splitbarcode-demux':
                splitbarcode = True
        if typearg and not typearg in ['illumina', 'mgig400', 'mgit7']:
            _raise(SanitycheckError("Unsupported protocol '" + typearg + "'"))
        if flag and not typearg == 'mgit7':
            log.info("Ignoring -f/--flag option because useless without '-t/--type mgit7'...")
        if typearg == 'illumina' and splitbarcode:
            log.info("Ignoring --splitbarcode-demux option because useless with '-t/--type illunina'...")

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
    def run_id(self):
        """
        For Ilumina runs :
         The run id from the run folder.
         Supports both default folder name configuration and GQ's globaly unique name convention.
        For MGI runs :
         The RUN ID from the run folder or from parameter
        """
        if self.args.type == 'illumina':
            if not hasattr(self, "_run_id"):
                if re.search(".*_\d+HS\d\d[AB]", self.run_dir):
                    m = re.search(".*/(\d+_[^_]+_\d+_[^_]+_(\d+)HS.+)", self.run_dir)
                    self._run_id = m.group(2)
                elif re.search(".*\d+_[^_]+_\d+_.+", self.run_dir):
                    m = re.search(".*/(\d+_([^_]+_\d+)_.*)", self.run_dir)
                    self._run_id = m.group(2)
                else:
                    log.warn("Unsupported folder name: " + self.run_dir)
        else:
            if not hasattr(self, "_run_id"):
                if self.args.run_id:
                    self._run_id = self.args.run_id
                else:
                    rundir_basename = os.path.basename(self.run_dir.rstrip('/'))
                    if (self.args.type == 'mgig400') and ("_" in rundir_basename):
                        [junk_food, self._run_id] = rundir_basename.split("_")
                    elif self.args.type == 'mgit7':
                        self._run_id = rundir_basename
                    else:
                        _raise(SanitycheckError("Error: Run ID could not be parsed from the RUN folder : " + self.run_dir))
        return self._run_id

    @property
    def run_dir(self):
        """
        The RUN ID from the run folder or from parameter
        """
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
        """
        The flow cell ID from the readset file
        """
        if not hasattr(self, "_flowcell_id"):
            all_flowcells = [readset.flow_cell for lane in self.lanes for readset in self.readsets[lane]]
            flowcells = set(all_flowcells)
            if len(flowcells) == 1:
                self._flowcell_id = flowcells.pop()
            elif len(flowcells) == 0:
                _raise(SanitycheckError("Error: No flowcell ids could be found."))
            else:
                # More than one flowcell found in the run, report which flowcells were found in which lanes.
                flowcell_counts = Counter(all_flowcells)
                readset_lanes = [lane for lane in self.lanes for _ in self.readsets[lane]]
                flowcell_found_in_lanes = {unique_flowcell: {lane for (flowcell, lane) in zip(all_flowcells, readset_lanes) if flowcell == unique_flowcell} for unique_flowcell in flowcell_counts}
                msg = ", ".join(["Flowcell '{}' in lanes: {}".format(flowcell, ",".join(sorted(lanes))) for (flowcell, lanes) in flowcell_found_in_lanes.items()])
                _raise(SanitycheckError("Error: Multiple flowcells found in this run. Please check the readset file. " + msg))
        return self._flowcell_id

    @property
    def raw_flag_dir(self):
        if self.args.raw_flag_dir:
            return self.args.raw_flag_dir
        else:
            _raise(SanitycheckError("Error: missing '-f/--flag' option!"))

    @property
    def raw_fastq_prefix(self):
        """
        """
        if not hasattr(self, "_raw_fastq_prefix"):
            if self.flowcell_id == self.run_id:
                self._raw_fastq_prefix = self.flowcell_id
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
            log.debug("yaaaaaa")
            self._json_flag_files = {}
            for lane in self.lanes:
                json_flag_file = os.path.join(self.output_dir, "flag." + lane + ".json")
                for filename in os.listdir(self.raw_flag_dir):
                    if re.match(self.run_id + "_" + lane + "_.+json", filename):
                        if not os.path.exists(json_flag_file):
                            if not os.path.exists(os.path.dirname(json_flag_file)):
                                os.makedirs(os.path.dirname(json_flag_file))
                            shutil.copy(os.path.join(self.raw_flag_dir, filename), json_flag_file)
                        self._json_flag_files[lane] = json_flag_file
                        log.info("JSON FLAG file for lane " + lane + " : " + json_flag_file)
                        break
                else:
                    _raise(SanitycheckError("Could not find any proper JSON flag file in " + self.raw_flag_dir + " for RUN " + self.run_id))
        log.debug("yooooo")
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

    @property
    def no_index_fastq(self):
        if not hasattr(self, "_no_index_fastq"):
            self._no_index_fastq = False
            if not self.args.type == 'illumina' and self.args.splitbarcode_demux:
                self._no_index_fastq = True
        return self._no_index_fastq

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
                if self.read2cycles[lane] == '0':
                    self._is_paired_end[lane] = False
                elif self.read2cycles[lane]:
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
                log.debug(self.get_index1cycles(lane))
        return self._index1cycles

    @property
    def index2cycles(self):
        if not hasattr(self, "_index2cycles"):
            self._index2cycles = {}
            for lane in self.lanes:
                self._index2cycles[lane] = self.get_index2cycles(lane)
        return self._index2cycles

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
    def index_per_readset(self):
        if not hasattr(self, "_index_per_readset"):
            return ""
        # Define in generate_clarity_sample_sheet()
        return self._index_per_readset

    @property
    def read_infos(self):
        if not hasattr(self, "_read_infos"):
            self._read_infos = self.parse_run_info_file()
        return self._read_infos

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
                if config.param('fastq_illumina', 'merge_undetermined', required=False, param_type='boolean'):
                    self._merge_undetermined[lane] = config.param('fastq_illumina', 'merge_undetermined')
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
                    "version" : "3.0",
                    "run" : self.run_id,
                    "instrument" : self.instrument,
                    "flowcell" : self.flowcell_id,
                    "lane" : lane,
                    "seqtype" : self.seqtype,
                    "sequencing_method" : "PAIRED_END" if self.is_paired_end[lane] else "SINGLE_END",
                    "steps" : [],
                    "readsets" : dict([
                        (
                            readset.name,
                            {
                                "sample_name": readset.sample.name,
                                "barcodes": readset.indexes,
                                "species": readset.species,
                                "reported_sex": readset.gender,
                                "fastq_1": readset.fastq1,
                                "fastq_2": readset.fastq2 if self.is_paired_end[lane] else None,
                                "bam": readset.bam + ".bam" if readset.bam else None,
                                "bai": readset.bam + ".bai" if readset.bam else None
                            }
                        ) for readset in self.readsets[lane]
                    ])
                }
        return self._report_hash

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

    def basecall(self):
        """
        Use write_fastq software from MGI to perform the base calling.
        Takes the raw .cal files from the sequencer and produces fastq files.
        Demultiplexing with MGI splitBarcode while doing the basecalling can
        be perform if requested with --splitbarcode-demux
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
                # Add the barcodes in the JSON flag file
                self.edit_mgi_t7_flag_file(lane)

                basecall_outputs, postprocessing_jobs = self.generate_basecall_outputs(lane)
                basecall_outputs.extend(
                    [
                        os.path.join(basecall_dir, self.run_id, "L0" + lane),
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + ".summaryReport.html"),
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + ".heatmapReport.html"),
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, "summaryTable.csv"),
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, "SequenceStat.txt"),
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, "BarcodeStat.txt")
                    ]
                )

                lane_jobs.append(
                    concat_jobs(
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
                            )
                        ],
                        name="basecall." + self.run_id + "." + lane,
                        samples=self.samples[lane],
                        report_files=[os.path.join(basecall_dir, self.run_id, "L0" + lane, "SequenceStat.txt")]
                    )
                )
                for readset in self.readsets[lane]:
                    readset.report_files['basecall'] = [os.path.join(basecall_dir, self.run_id, "L0" + lane, "SequenceStat.txt")]

                if postprocessing_jobs:
                    jobs_to_throttle.extend(postprocessing_jobs)

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
                if int(self.index1cycles[lane]) + int(self.index2cycles[lane]) == 0 and len(self.readsets[lane]) > 1:
                    err_msg = "LANE SETTING ERROR :\n"
                    err_msg += "Unable to demultiplex " + str(len(self.readsets[lane])) + " samples : No barcode in fastq files...\n(in "
                    err_msg += self.run_dir + ")"
                    _raise(SanitycheckError(err_msg))

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
                    concat_jobs(
                        [
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
                        ],
                        name="bcl2fastq_index." + self.run_id + "." + lane,
                        samples=self.samples[lane],
                        removable_files=[os.path.join(self.output_dir, "index", "L00" + lane)],
                        report_files=[os.path.join(self.output_dir, "index", "L00" + lane, "Stats", "Stats.json")]
                    )
                )
                for readset in self.readsets[lane]:
                    readset.report_files['index'] = [os.path.join(self.output_dir, "index", "L00" + lane, "Stats", "Stats.json")]

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

                self.add_to_report_hash("index",lane,  lane_jobs)
                self.add_copy_job_inputs(lane_jobs, lane)
                jobs.extend(lane_jobs)

        return jobs

    def fastq(self):
        """
        For Illumina

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

        For MGI-G400

            First copy all the files of the lane from the sequencer deposit folder
            to the processing folder, into "raw_fastq".
            Then, perform demultplexing of the reads with fgbio DemuxFastqs

        For MGI-T7

            Perform demultiplexing of the reads with fgbio DemuxFastqs
            (skipped with --splitbarcode-demux)
        """
        if self.args.type == 'illumina':
            jobs = self.fastq_illumina()
        else:
            jobs = self.mgi_fastq()
        return jobs

    def fastq_illumina(self):
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
            fastq_outputs, final_fastq_jobs = self.generate_bcl2fastq_outputs(lane)
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
                            mask=self.mask[lane],
                            ini_section='fastq_illumina'
                        ),
                        run_processing_tools.bcl2fastq(
                            input,
                            fastq_outputs,
                            output_dir_noindex,
                            casava_sample_sheet_noindex,
                            self.run_dir,
                            lane,
                            self.bcl2fastq_extra_option,
                            ini_section='fastq_illumina'
                    )],
                    name="fastq_illumina." + self.run_id + "." + lane,
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
                    mask=self.mask[lane],
                    ini_section='fastq_illumina'
                )
                bcl2fastq_job.name = "fastq_illumina." + self.run_id + "." + lane
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
#                    number_of_mismatches=self.number_of_mismatches,
                    lane_number=self.lane_number,
#                    mask=self.mask[lane],
#                    technology=config.param('fastq', 'technology'),
#                    run_id=self.run_id
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
#                    technology=config.param('fastq', 'technology'),
#                    run_id=self.run_id
                )
                job = Job(
                    fastq_outputs,
                    ["notificationFastqEnd." + lane + ".out"],
                    command=notification_command_end,
                    name="fastq_notification_end." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
#                lane_jobs.append(job)

            self.add_to_report_hash("fastq", lane, lane_jobs)
            jobs.extend(lane_jobs)

        return jobs


    def mgi_fastq(self):
        """
        *** In the future, may generate the fastq files from the raw CAL files. ***
        Perform demultplexing of the reads with fgbio DemuxFastqs
        """

        ini_section = None
        if self.args.type == 'mgig400':
            ini_section = 'fastq_g400'
        elif self.args.type == 'mgit7':
            ini_section = 'fastq_t7'
        if not ini_section:
            _raise(SanitycheckError("Could not determine which section to use for fastq step from given protocol " + self.protocol))
        jobs = []

        if self.args.splitbarcode_demux:
            log.info("Demultiplexing done during the basecalling... Skipping fastq step...")

        else:
            log.info("No demultiplexing done yet... Processing fastq step...")

            for lane in self.lanes:

                if int(self.index1cycles[lane]) + int(self.index2cycles[lane]) == 0 and len(self.readsets[lane]) > 1:
                     err_msg = "LANE SETTING ERROR :\n"
                     err_msg += "Unable to demultiplex " + str(len(self.readsets[lane])) + " samples : No barcode in fastq files...\n(in "
                     err_msg += self.run_dir + ")"
                     _raise(SanitycheckError(err_msg))

                elif int(self._index1cycles[lane]) + int(self._index2cycles[lane]) == 0:
                    log.info("No barcode cycles in the lane... Skipping fastq step for lane " + lane + "...")

                else:

                    lane_jobs = []

                    if ini_section == 'fastq_g400':
                        # For MGI G400 runs, raw reads are in non demultiplexed fastq files
                        # Here is all the prep to copy the lane folder from the sequencer deposit folder into raw_fastq_dir
                        unaligned_dir = os.path.join(self.output_dir, "Unaligned." + lane)
                        raw_fastq_dir = os.path.join(unaligned_dir, "raw_fastq")
                        copy_done_file = os.path.join(raw_fastq_dir, "copy_done.Success")

                        raw_name_prefix = self.raw_fastq_prefix +  "_L0" + lane

                        # Start setting the copy jobs output files
                        copy_job_output_dependency = [
                            os.path.join(raw_fastq_dir, raw_name_prefix + ".summaryReport.html"),
                            os.path.join(raw_fastq_dir, raw_name_prefix + ".heatmapReport.html"),
                            os.path.join(raw_fastq_dir, "summaryTable.csv")
                        ]

                        if self.readsets[lane][0].run_type == "PAIRED_END":
                            copy_job_output_dependency.append(os.path.join(raw_fastq_dir, raw_name_prefix + "_read_1.fq.gz"))
                            copy_job_output_dependency.append(os.path.join(raw_fastq_dir, raw_name_prefix + "_read_2.fq.gz"))
                            copy_job_output_dependency.append(os.path.join(raw_fastq_dir, raw_name_prefix + "_read_1.fq.fqStat.txt"))
                            copy_job_output_dependency.append(os.path.join(raw_fastq_dir, raw_name_prefix + "_read_2.fq.fqStat.txt"))
                        else:
                            copy_job_output_dependency.append(os.path.join(raw_fastq_dir, raw_name_prefix + "_read.fq.gz"))
                            copy_job_output_dependency.append(os.path.join(raw_fastq_dir, raw_name_prefix + "_read.fq.fqStat.txt"))

                        copy_job_output_dependency.append(os.path.join(raw_fastq_dir, raw_name_prefix + "_read.report.html"))

                        # If needed, set a renaming job here
                        rename_job = None
                        if (len(self.readsets[lane]) == 1) or (int(self._index1cycles[lane]) + int(self._index2cycles[lane]) == 0):
                            readset = self.readsets[lane][0]
                            if readset.run_type == "PAIRED_END":
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
                                        ),
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
                            else:
                                rename_job = concat_jobs(
                                    [
                                        bash.mkdir(os.path.dirname(readset.fastq1)),
                                        bash.cp(
                                            os.path.join(raw_fastq_dir, raw_name_prefix + "_read.fq.gz"),
                                            readset.fastq1,
                                        ),
                                        bash.cp(
                                            os.path.join(raw_fastq_dir, raw_name_prefix + "_read.fq.fqStat.txt"),
                                            re.sub("", "fqStat.txt", readset.fastq1),
                                        ),
                                        bash.cp(
                                            os.path.join(raw_fastq_dir, raw_name_prefix + "_read.report.html"),
                                            re.sub("_R1_001.fastq.gz", "_read.report.html", readset.fastq1),
                                        )
                                    ]
                                )
                            rename_job.output_files.append(copy_done_file)
                            rename_job.name = ini_section + ".rename." + self.run_id + "." + lane
                            rename_job.samples = self.samples[lane]

                        # Set the copy job along with some md5 checksum (unless already done during previous execution of the pipeline...)
                        if not os.path.exists(copy_done_file) or self.force_jobs:
                            copy_job = concat_jobs(
                                [
                                    bash.mkdir(raw_fastq_dir),
                                    bash.cp(
                                        os.path.join(self.run_dir, "L0" + lane, "."),
                                        raw_fastq_dir,
                                        recursive=True
                                    ),
                                    pipe_jobs(
                                        [
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
                                        ]
                                    ),
                                    pipe_jobs(
                                        [
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
                                        ]
                                    ),
                                    bash.md5sum(
                                        copy_done_file + ".md5",
                                        None,
                                        check=True
                                    ),
                                    bash.touch(copy_done_file)
                                ],
                                input_dependency=[os.path.join(self.run_dir, "L0" + lane, ".")],
                                output_dependency=copy_job_output_dependency + [re.sub(os.path.dirname(self.bioinfo_files[lane]), raw_fastq_dir, self.bioinfo_files[lane]), copy_done_file],
                                name=ini_section + ".copy_raw." + self.run_id + "." + lane,
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

                    jobs_to_throttle = []

                    input_fastq_dir = os.path.join(self.output_dir, "Unaligned." + lane, "raw_fastq")

                    demuxfastqs_outputs, postprocessing_jobs = self.generate_demuxfastqs_outputs(lane)

                    tmp_output_dir = os.path.dirname(demuxfastqs_outputs[0])
                    tmp_metrics_file = os.path.join(tmp_output_dir, self.run_id + "." + lane + ".DemuxFastqs.metrics.txt")
                    metrics_file = os.path.join(self.output_dir, "Unaligned." + lane, self.run_id + "." + lane + ".DemuxFastqs.metrics.txt")
                    demuxfastqs_outputs.append(metrics_file)

                    if ini_section == 'fastq_g400':
                        raw_name_prefix = self.raw_fastq_prefix +  "_L0" + lane
                        input1 = os.path.join(input_fastq_dir, raw_name_prefix + "_read_1.fq.gz")
                        input2 = os.path.join(input_fastq_dir, raw_name_prefix + "_read_2.fq.gz")
                        if self.readsets[lane][0].run_type == "PAIRED_END":
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
                        else:
                            input1 = os.path.join(input_fastq_dir, raw_name_prefix + "_read.fq.gz")
                            demultiplex_job = run_processing_tools.demux_fastqs_single_end(
                                os.path.join(self.output_dir, "samplesheet." + lane + ".csv"),
                                self.number_of_mismatches,
                                self.mask[lane],
                                demuxfastqs_outputs,
                                tmp_metrics_file,
                                input1,
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
                            output_dependency=demuxfastqs_outputs,
                            report_files=[metrics_file]
                        )
                    )
                    for readset in self.readsets[lane]:
                        readset.report_files['fastq'] = [metrics_file]

                    if postprocessing_jobs:
                        jobs_to_throttle.extend(postprocessing_jobs)

                    if self.args.type == 'mgit7':
                        lane_jobs.extend(jobs_to_throttle)
                    else:
                        lane_jobs.extend(self.throttle_jobs(jobs_to_throttle, lane))

                    self.add_to_report_hash("fastq", lane, lane_jobs)
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
                    concat_jobs(
                        [
                            bash.mkdir(output_dir),
                            bvatools.readsqc(
                                readset.fastq1,
                                readset.fastq2,
                                "FASTQ",
                                region_name,
                                output_dir
                            )
                        ],
                        name="qc_graphs." + readset.name + ".qc." + self.run_id + "." + lane,
                        samples=[readset.sample],
                        report_files=[os.path.join(output_dir, "mpsQC_" + region_name + "_stats.xml")]
                    )
                )
                readset.report_files['qc_graphs'] = [os.path.join(output_dir, "mpsQC_" + region_name + "_stats.xml")]

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

        if self.args.type == 'mgit7':
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
                    report_file = None
                    if input1 == readset.fastq1:
                        job_suffix = "R1."
                        input2 = None
                        unzip = True
                        if not len(self.readsets[lane]) == 1:
                            report_file = outputs[2]
                            readset.report_files['fastqc'] = [report_file]
                    elif input1 == readset.fastq2:
                        job_suffix = "R2."
                        input2 = None
                        unzip=True
                    elif input1 == readset.index_fastq1:
                        job_suffix = "Barcodes."
                        input2 = readset.index_fastq2 if self.is_dual_index[lane] else None
                    lane_jobs.append(
                        concat_jobs(
                            [
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
                                )
                            ],
                            name="fastqc." + readset.name + "_" + job_suffix + "." + self.run_id + "." + lane,
                            samples=[readset.sample],
                            report_files=[report_file]
                        )
                    )

            self.add_to_report_hash("fastqc", lane, lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)

        if self.args.type == 'mgit7':
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
                    concat_jobs(
                        [
                            bash.mkdir(os.path.dirname(output_json_path), remove=True),
                            fastp.fastp_basic_qc(readset.fastq1, readset.fastq2, output_json_path, output_html_path),
                        ],
                        name="fastp." + readset.name + "." + self.run_id + "." + lane,
                        samples=[readset.sample],
                        report_files=[output_json_path]
                    )
                )
                readset.report_files['fastp'] = [output_json_path]
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

        nb_blast_to_do = config.param('blast', 'nb_blast_to_do', param_type="posint")
        is_nb_blast_per_lane = config.param('blast', 'is_nb_for_whole_lane', param_type="boolean")

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
                    samples=[readset.sample],
                    report_files=[output]
                )
                lane_jobs.append(job)
                readset.report_files['blast'] = [output]

            self.add_to_report_hash("blast", lane, lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)

        if self.args.type == 'mgit7':
            return jobs
        else:
            return self.throttle_jobs(jobs)

    # Not used anymore...
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
                    samples=[readset.sample],
                    report_files=[os.path.join(ouput_directory, readset.name + ".sample_tag_stats.csv")]
                ))
                readset.report_files['sample_tag'] = [os.path.join(ouput_directory, readset.name + ".sample_tag_stats.csv")]

            self.add_to_report_hash("sample_tag", lane, lane_jobs)
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

        if self.args.type == 'mgit7':
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

                job = picard.mark_duplicates(
                    [input],
                    output,
                    metrics_file
                )
                job.name = "picard_mark_duplicates." + readset.name + ".dup." + self.run_id + "." + lane
                job.samples = [readset.sample]
                job.report_files = [metrics_file]
                lane_jobs.append(job)
                readset.report_files['picard_mark_duplicates'] = [metrics_file]

            self.add_to_report_hash("picard_mark_duplicates", lane, lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)

            jobs.extend(lane_jobs)

        if self.args.type == 'mgit7':
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
                report_files = []
                for job in job_list:
                    job.samples = [readset.sample]
                    report_files.extend(job.report_files)
                lane_jobs.extend(job_list)
                readset.report_files['metrics'] = report_files

            self.add_to_report_hash("metrics", lane, lane_jobs)
            self.add_copy_job_inputs(lane_jobs, lane)
            jobs.extend(lane_jobs)

        if self.args.type == 'mgit7':
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

            if config.param('md5', 'one_job', required=False, param_type="boolean"):
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
        Generate a JSON file reporting the whole pipeline.
        The jobs of this step actually update the JSON report as the pipeline is running
        """
        jobs = []

        for lane in self.lanes:
            lane_jobs = []

            report_dir = os.path.join(self.output_dir, "report")

            self.generate_lane_json_report_file(lane)

            # metrics to JSON
            # Loop over all the steps of the pipeline
            step_list = [step for step in self.step_list if step.jobs]
            for step in step_list:
                report_step_jobs = []
                if step.name in ['basecall', 'fastq', 'index']:
                    step_report_files = list(set([report_file for readset in self.readsets[lane] for report_file in readset.report_files[step.name] if step.name in readset.report_files]))
                    if step_report_files:
                        report_job = tools.run_processing_metrics_to_json(
                             self.run_validation_report_json[lane],
                             step.name,
                             self.args.type,
                             step_report_files
                        )
                        report_job.name = "report." + step.name + "." + self.run_id + "." + lane
                        report_job.samples = self.samples[lane]
                        report_step_jobs.append(report_job)
                else:
                    for readset in self.readsets[lane]:
                        if step.name in readset.report_files:
                            report_job = tools.run_processing_metrics_to_json(
                                self.run_validation_report_json[lane],
                                step.name,
                                self.args.type,
                                readset.report_files[step.name],
                                readset.name
                            )
                            report_job.name = "report." + step.name + "." + readset.name + "." + self.run_id + "." + lane
                            report_job.samples = self.samples[lane]
                            report_step_jobs.append(report_job)

                lane_jobs.extend(self.throttle_jobs(report_step_jobs, f"{step.name}.{lane}"))

                # checkpoint file
                step_checkpoint_file = os.path.join(self.job_output_dir, "checkpoint", step.name + "." + self.run_id + "." + lane + ".stepDone")
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
        Copy the whole processing folder to where they can be serve or loaded into a LIMS
        """
        jobs_to_concat = []

        if not self.args.type == 'illumina':
            full_destination_folder = os.path.join(
                config.param("copy", "destination_folder", param_type="dirpath"),
                self.seq_category,
                self.year,
                self.date + "_" + self.instrument + "_" + self.run_number + "_" + self.flowcell_position + self.flowcell_id + "_" + self.sequencer_run_id + "-" + self.seqtype
            )
        else:
            full_destination_folder = os.path.join(
                config.param("copy", "destination_folder", param_type="dirpath"),
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

            exclude_bam = config.param('copy', 'exclude_bam', required=False, param_type='boolean')
            exclude_fastq_with_bam = config.param('copy', 'exclude_fastq_with_bam', required=False, param_type='boolean')
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

    # Obsolete...
    def end_copy_notification(self):
        """
        Send an optional notification to notify that the copy is finished.

        The command used is in the configuration file. This step is skipped when no
        command is provided.
        """
        jobs = []

        full_destination_folder = os.path.join(
            config.param('copy', 'destination_folder', param_type="dirpath"),
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

    def final_notification(self):
        """
        Writes a simple '.done' file when the whole pipeline is over
        """
        jobs = []

        inputs = []
        for lane in self.lanes:
            if not self.args.type == 'illumina':
                full_destination_folder = os.path.join(
                    config.param("copy", "destination_folder", param_type="dirpath"),
                    self.seq_category,
                    self.year,
                    self.date + "_" + self.instrument + "_" + self.run_number + "_" + self.flowcell_position + self.flowcell_id + "_" + self.sequencer_run_id + "-" + self.seqtype
                )
            else:
                full_destination_folder = os.path.join(
                    config.param("copy", "destination_folder", param_type="dirpath"),
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

    def get_read1cycles(self, lane):
        """
        Parse the BioInfo.csv file for the number of read 1 cycles in the run
        """
        if self.args.type == "illumina":
            [read1cycles, read2cycles] = self.get_illumina_readcycles()
            return read1cycles
        elif self.args.type == "mgig400":
            return self.bioinfo_hash[lane]["Read1 Cycles"]

        elif self.args.type == "mgit7":
            return self.json_flag_hash[lane]['Read1']

        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_read2cycles(self, lane):
        """
        Parse the BioInfo.csv file for the number of read 2 cycles in the run
        """
        if self.args.type == "illumina":
            [read1cycles, read2cycles] = self.get_illumina_readcycles()
            return read2cycles
        elif self.args.type == "mgig400":
            return self.bioinfo_hash[lane]["Read2 Cycles"]

        elif self.args.type == "mgit7":
            return self.json_flag_hash[lane]['Read2']
        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_illumina_readcycles(self):
        """
        Returns the number of cycles for each reads of the run.
        """
        for read in [read for read in self.read_infos if not read.is_index]:
            if read.number in [1, 2]:
                read1cycles = read.nb_cycles
            if read.number in [3, 4]:
                read2cycles = read.nb_cycles
        return [read1cycles, read2cycles]

    def get_index1cycles(self, lane):
        """
        Parse the BioInfo.csv file for the number of index 1 cycles in the run
        """
        if self.args.type == "illumina":
            [index1cycles, index2cycles] = self.get_illumina_indexcycles()
            return index1cycles
        elif self.args.type == "mgig400":
            if self.bioinfo_hash[lane]["Barcode"] == '':
                return "0"
            else:
                return self.bioinfo_hash[lane]["Barcode"]
            _raise(SanitycheckError("Could not get Index 1 Cycles from " + self.bioinfo_files[lane]))

        elif self.args.type == "mgit7":
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
        if self.args.type == "illumina":
            [index1cycles, index2cycles] = self.get_illumina_indexcycles()
            return index2cycles
        elif self.args.type == "mgig400":
            if self.bioinfo_hash[lane]["Dual Barcode"] == '':
                return "0"
            else:
                return self.bioinfo_hash[lane]["Dual Barcode"]

        elif self.args.type == "mgit7":
            if self.json_flag_hash[lane]["Dual Barcode"] == '':
                return "0"
            else:
                return self.json_flag_hash[lane]["Dual Barcode"]
        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_illumina_indexcycles(self):
        """
        Returns the number of cycles for each reads of the run.
        """
        for read in [read for read in self.read_infos if read.is_index]:
            if read.number in [1, 2]:
                index1cycles = read.nb_cycles
            if read.number in [3, 4]:
                index2cycles = read.nb_cycles
        return [index1cycles, index2cycles]

    def get_sequencer_index_length(self):
        """
        Returns the total number of index cycles of the run.
        """
        return sum(index_read.nb_cycles for index_read in [read for read in self.read_infos if read.is_index])

    def get_sequencer_minimum_read_length(self):
        """
        Returns the minimum number of cycles of a real read (not indexed).
        """
        return min(read.nb_cycles for read in [read for read in self.read_infos if not read.is_index])

    # Obsolete...
    def has_single_index(self, lane):
        """
        Returns True when there is at least one sample on the lane that doesn't use double-indexing or we only have
        one read of indexes.
        """
        return len([readset for readset in self.readsets[lane] if ("-" not in readset.index)]) > 0 or\
               len([read for read in self.read_infos[lane] if read.is_index]) < 2

    def get_smallest_index_length(self, lane):
        if self.args.type == 'illumina':
            return self.get_illumina_smallest_index_length(lane)
        else:
            return self.get_mgi_smallest_index_length(lane)

    def get_illumina_smallest_index_length(self, lane):
        """
        Returns a list (for each index read of the lane) of the minimum between the number of index cycle on the
        sequencer and all the index lengths.
        """
        run_index_lengths = [r.nb_cycles for r in self.read_infos if r.is_index] # from RunInfo

#        if len(run_index_lengths) == 0 and len(self.readsets[lane]) > 1:
#            _raise(SanitycheckError("Multiple samples on lane '" + lane + "', but no indexes were read from the sequencer."))

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

    def get_mgi_smallest_index_length(self, lane):
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

    def get_instrument(self, lane):
        """
        For Ilumina runs :
         Parse the RunInfo.xml file for the instrument name the run has been running on
        For MGI G400 runs :
          Parse the BioInfo.csv file for the instrument name the run has been running on
        For MGI T7 runs :
          Parse the JSON flag file for the instrument name the run has been running on
        """
        if self.args.type == "illumina":
            return Xml.parse(os.path.join(self.run_dir, "RunInfo.xml")).getroot().find('Run').find('Instrument').text
        elif self.args.type == "mgig400":
            return self.bioinfo_hash[lane]["Machine ID"]

        elif self.args.type == "mgit7":
            return self.json_flag_hash[lane]['Machine ID']
        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_flowcell_position(self, lane):
        """
        For Ilumina runs :
          Not supposed to be called because not needed
        For MGI G400 runs :
          Parse the BioInfo.csv file for flowcell position of the run
        For MGI T7 runs :
          Parse the JSON flag file for flowcell position of the run
        """
        if self.args.type == "illumina":
            _raise(SanitycheckError("Wrong call of get_flowcell_position..."))
        elif self.args.type == "mgig400":
            return self.bioinfo_hash[lane]["Flowcell Pos"]

        elif self.args.type == "mgit7":
            return self.json_flag_hash[lane]["Flow Cell Pos"]

        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_sequencer_run_id(self, lane):
        """
        For Ilumina runs :
          Not supposed to be called because not needed
        For MGI G400 runs :
          Parse the BioInfo.csv file for the run id given by the sequencer
        For MGI T7 runs :
          Parse the JSON flag file for the run id given by the sequencer
        """
        if self.args.type == "illumina":
            _raise(SanitycheckError("Wrong call of get_sequencer_run_id..."))
        elif self.args.type == "mgig400":
            # dnb_id format looks like : 10074MG01B_Lane4
            # where run_id is : 10074MG01B
            return self.bioinfo_hash[lane]["DNB ID"].split("_")[0]

        elif self.args.type == "mgit7":
            return self.json_flag_hash[lane]['Flow Cell ID'].split("_")[0]

        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_run_number(self, lane):
        """
        For Ilumina runs :
          Not supposed to be called because not needed
        For MGI G400 runs :
          Parse the BioInfo.csv file for the run number
        For MGI T7 runs :
          Parse the JSON flag file for the run number
        """
        if self.args.type == "illumina":
            _raise(SanitycheckError("Wrong call of get_run_number..."))
        elif self.args.type == "mgig400":
            # dnb_id format looks like : 10074MG01B_Lane4
            # where run_id is : 10074MG01B
            # and run counter is : 10074
            run_number = self.bioinfo_hash[lane]["DNB ID"].split("_")[0][:-5]
            # sometimes, format is misleading : 1074
            # so we correct it to : 10074
            while len(run_number) < 5:
                run_number[:1] + "0" +  run_number[1:]
            return run_number

        elif self.args.type == "mgit7":
            return self.json_flag_hash[lane]['Flow Cell ID'].split("_")[0][-5:]

        else:
            _raise(SanitycheckError("Unknown protocol : " + self.args.type))

    def get_seqtype(self):
        """
        Determine which kind of sequencing (iseq, miseq, novaseq, hiseqx, hiseq4000, hiseq2500 or dnbseqg400) was performed,
        depending on the instrument used for the run
        """

        instrument = self.instrument
        instrument_file = config.param('DEFAULT', 'instrument_list_file', param_type='filepath', required=False)
        if not (instrument_file and os.path.isfile(instrument_file)):
            instrument_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'instrument_list.csv')

        return subprocess.check_output("grep -m1 '"+instrument+"' %s | awk -F',' '{print $3}'" % instrument_file, shell=True, text=True).strip()

    def get_seq_category(self):
        if "hiseq" in self.seqtype:
            return "hiseq"
        else:
            return self.seqtype

    def validate_barcodes(self, lane):
        if self.args.type == "illumina":
            self.illumina_validate_barcodes(lane)
        else:
            self.mgi_validate_barcodes(lane)

    def illumina_validate_barcodes(self, lane):
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

    def mgi_validate_barcodes(self, lane):
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
        if self.args.type == "illumina":
            return self.illumina_get_mask(lane)
        else:
            return self.mgi_get_mask(lane)

    def illumina_get_mask(self, lane):
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

    def mgi_get_mask(self, lane):
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

        mask = self.read1cycles[lane] + 'T'
        if self.read2cycles[lane] != '0':
            mask += ' ' + self.read2cycles[lane] + 'T'

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

    def edit_mgi_t7_flag_file(self, lane):
        """
        insert the barcode information in the flag file output from the sequencer
        This Flag file is used to call splitBarcode
        """

        json_flag_file = self.json_flag_files[lane]

        # index_lengths = self.get_smallest_index_length(lane)

        # get the barcode names & sequences to add in the JSON flag file
        all_indexes = {}
        for readset in self.readsets[lane]:
            # for index_dict in readset.index:
            #     all_indexes[readset.index_name] = index_dict
            for readset_index in readset.indexes:
                all_indexes[readset.index_name] = readset_index
        with open(json_flag_file, 'r') as json_fh:
            json_flag_content = json.load(json_fh)
        if self.is_dual_index[lane]:
            json_flag_content['speciesBarcodes'] = dict([(index_name, str(Seq(readset_index['INDEX2'] + readset_index['INDEX1']).reverse_complement())) for index_name, readset_index in all_indexes.items()])
        else:
            json_flag_content['speciesBarcodes'] = dict([(index_name, str(Seq(readset_index['INDEX1']).reverse_complement())) for index_name, readset_index in all_indexes.items()])
        json_flag_content['barcodeStartPos'] = json_flag_content['TotalCycle'] - json_flag_content['barcodeLength'] + 1
        json_flag_content['SpeciesMismatch'] = self.number_of_mismatches

        with open(json_flag_file, 'w') as out_json_fh:
            json.dump(json_flag_content, out_json_fh, indent=4)

        log.info("BARCODES added in FLAG file : " + json_flag_file)

    def generate_lane_sample_sheet(self, lane):
        if self.args.type == 'illumina':
            self.generate_illumina_lane_sample_sheet(lane)
        else:
            self.generate_mgi_lane_sample_sheet(lane)

    def generate_illumina_lane_sample_sheet(self, lane):
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
            open(csv_file, 'w'),
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
        if re.search("I", mask) and not self.args.allow_barcode_collision:
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
                _raise(SanitycheckError("HaloPlex libraries cannot be mixed with DUAL INDEX libraries"))
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
        final_index1 = self.index1cycles[lane] if overindex1 is None else overindex1
        final_index2 = self.index2cycles[lane] if overindex2 is None else overindex2

        self._mask[lane] = config.param('fastq_illumina', 'overmask') if config.param('fastq_illumina', 'overmask', required=False, param_type='string') else final_mask
        self._index1cycles[lane] = config.param('fastq_illumina', 'overindex1') if config.param('fastq_illumina', 'overindex1', required=False, param_type='int') else final_index1
        self._index2cycles[lane] = config.param('fastq_illumina', 'overindex2') if config.param('fastq_illumina', 'overindex2', required=False, param_type='int') else final_index2

        # If the second index exists
        if self.index2cycles != 0:
            dualindex_demultiplexing = True

        # In case of HaloPlex-like masks, R2 is actually I2 (while R3 is R2),
        # Proceed to dualindex demultiplexing
        if ''.join(i for i in self.mask[lane] if not i.isdigit()) == "Y,I,Y,Y":
            dualindex_demultiplexing = True
        output_dir = os.path.join(self.output_dir, "Unaligned." + lane)
        casava_sample_sheet = os.path.join(self.output_dir, "casavasheet." + lane + ".indexed.csv")

        count = 0
        index_per_readset = {}
        index_lengths = self.get_smallest_index_length(lane)
        for readset in self.readsets[lane]:
            count += 1
            index_per_readset[readset.name] = readset.indexes

            for idx, readset_index in enumerate(readset.indexes):
                # Barcode sequence should only match with the barcode cycles defined in the mask
                # so we adjust the lenght of the index sequences accordingly for the "Sample_Barcode" field
                if self.is_dual_index[lane]:
                    sample_barcode = readset_index['INDEX1'][0:index_lengths[0]] + readset_index['INDEX2'][0:index_lengths[1]]
                else:
                    sample_barcode = readset_index['INDEX1'][0:index_lengths[0]]
                if self.last_index < len(sample_barcode):
                    sample_barcode = sample_barcode[0:self.last_index]
                if self.first_index > 1:
                    sample_barcode = sample_barcode[self.first_index-1:]

                readset_index['BARCODE_SEQUENCE'] = sample_barcode
                readset.indexes[idx] = readset_index

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
                open(casava_sample_sheet_noindex, 'w'),
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
                        "pct_on_index_in_lane": None,
                        "pct_of_the_lane": None,
                        "pct_perfect_barcode": None,
                        "pct_one_mismatch_barcode": None,
                        "pf_clusters": None,
                        "yield": None,
                        "mean_quality_score": None,
                        "pct_q30_bases": None
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
                        "pf_read_alignment_rate": None,
                        "freemix": None,
                        "inferred_sex": None,
                        "adapter_dimers": None,
                        "mean_coverage": None,
                        "aligned_dup_rate": None,
                        "sex_concordance": None
                    }
                }
            )

        # Adding the MultiQC input file list
        self.report_hash[lane]["multiqc_inputs"] = []
        step_list = [step for step in self.step_list if step.jobs]
        for step in step_list:
            if step.name in ['basecall', 'fastq', 'index']:
                step_report_files = list(set([report_file for job in step.jobs for report_file in job.report_files for sample in job.samples for readset in sample.readsets if job.report_files and readset.lane == lane]))
                self.report_hash[lane]["multiqc_inputs"].extend([os.path.relpath(path, self.report_dir[lane]) for path in step_report_files])
            else:
                for readset in self.readsets[lane]:
                    step_report_files = list(set([report_file for job in step.jobs for sample in job.samples for sreadset in sample.readsets for report_file in job.report_files if job.report_files and sreadset.lane == lane and sreadset.name == readset.name]))
                    self.report_hash[lane]["multiqc_inputs"].extend([os.path.relpath(path, self.report_dir[lane]) for path in step_report_files])

        if not os.path.exists(os.path.dirname(self.run_validation_report_json[lane])):
             os.makedirs(os.path.dirname(self.run_validation_report_json[lane]))
        if not os.path.exists(self.run_validation_report_json[lane]):
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
                readset_r1_outputs.extend([
                    os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_" + index['INDEX_NAME'] + "_1.fq.gz"),
                    os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_" + index['INDEX_NAME'] + "_1.fq.fqStat.txt")
                ])
                if readset.run_type == "PAIRED_END":
                    readset_r2_outputs.extend([
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_" + index['INDEX_NAME'] + "_2.fq.gz"),
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_" + index['INDEX_NAME'] + "_2.fq.fqStat.txt")
                    ])

            # If True, then merge the 'Undetermined' reads
            if self.merge_undetermined[lane]:
                readset_r1_outputs.extend([
                    os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_undecoded_1.fq.gz"),
                    os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_undecoded_1.fq.fqStat.txt")
                ])
                if readset.run_type == "PAIRED_END":
                    readset_r2_outputs.extend([
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_undecoded_2.fq.gz"),
                        os.path.join(basecall_dir, self.run_id, "L0" + lane, self.raw_fastq_prefix +  "_L0" + lane + "_undecoded_2.fq.fqStat.txt")
                    ])
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
                                    self.awk_read_1_processing_command(lane, readset)
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
                                        self.awk_read_2_processing_command(lane, readset)
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
                        self.awk_read_1_processing_command(lane)
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
                            self.awk_read_2_processing_command(lane)
                        )
                    ],
                    output_dependency=outputs,
                    name="fastq_convert.R2.unmatched." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
            )

        return basecall_outputs, postprocessing_jobs

    def generate_bcl2fastq_outputs(self, lane):
        bcl2fastq_outputs = []
        final_fastq_jobs = []
        count = 0
        merge = self.merge_undetermined[lane]

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

    def generate_demuxfastqs_outputs(self, lane):
        demuxfastqs_outputs = []
        postprocessing_jobs = []
        output_dir = os.path.join(self.output_dir, "Unaligned." + lane)

        for readset in self.readsets[lane]:
            readset_r1_outputs = []
            readset_r2_outputs = []

            for index in readset.indexes:
                if readset.run_type == "PAIRED_END":
                    readset_r1_outputs.append(
                        os.path.join(output_dir, "tmp", index['SAMPLESHEET_NAME']+"-"+index['SAMPLESHEET_NAME']+'_'+index['INDEX_NAME']+"-"+index['BARCODE_SEQUENCE']+"_R1.fastq.gz")
                    )
                    readset_r2_outputs.append(
                        os.path.join(output_dir, "tmp", index['SAMPLESHEET_NAME']+"-"+index['SAMPLESHEET_NAME']+'_'+index['INDEX_NAME']+"-"+index['BARCODE_SEQUENCE']+"_R2.fastq.gz")
                    )
                else:
                    readset_r1_outputs.append(
                        os.path.join(output_dir, "tmp", index['SAMPLESHEET_NAME']+"-"+index['SAMPLESHEET_NAME']+'_'+index['INDEX_NAME']+"-"+index['BARCODE_SEQUENCE']+".fastq.gz")
                    )

            # If True, then merge the 'Undetermined' reads
            if self.merge_undetermined[lane]:
                if readset.run_type == "PAIRED_END":
                    readset_r1_outputs.append(
                        os.path.join(output_dir, "tmp", "unmatched_R1.fastq.gz")
                    )
                    readset_r2_outputs.append(
                        os.path.join(output_dir, "tmp", "unmatched_R2.fastq.gz")
                    )
                else:
                    readset_r1_outputs.append(
                        os.path.join(output_dir, "tmp", "unmatched.fastq.gz")
                    )


            # Processing R1 fastq outputs :
            #   convert headers from MGI to Illumina format using zcat and awk
            demuxfastqs_outputs.extend(readset_r1_outputs)
            outputs = [readset.fastq1]
            if not self.is_paired_end[lane]:
                outputs.append(readset.index_fastq1)
                if self.is_dual_index[lane]:
                     outputs.append(readset.index_fastq2)
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
                                    self.awk_read_1_processing_command(lane, readset)
                                ),
                                bash.gzip(
                                    None,
                                    readset.fastq1
                                ) if self.is_paired_end[lane] else None
                            ]
                        )
                    ],
                    output_dependency=outputs,
                    name="fastq_convert.R1." + readset.name + "." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
            )

            if readset.run_type == "PAIRED_END":
                # Processing "raw R2" fastq (also contains barcode sequences) :
                #   convert headers from MGI to Illumina format
                #   while extracting I1 and I2 to build clean R2, I1 and R2 fastq
                #   using zcat and awk
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
                                        self.awk_read_2_processing_command(lane, readset)
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
        if readset.run_type == "PAIRED_END":
            unmatched_R1_fastq = os.path.join(output_dir, "tmp", "unmatched_R1.fastq.gz")
        else:
            unmatched_R1_fastq = os.path.join(output_dir, "tmp", "unmatched.fastq.gz")

        if unmatched_R1_fastq not in demuxfastqs_outputs:
            demuxfastqs_outputs.append(unmatched_R1_fastq)
        outputs = [os.path.join(output_dir, "Undetermined_S0_L00" + lane + "_R1_001.fastq.gz")]
        if not self.is_paired_end[lane]:
            unaligned_i1 = os.path.join(output_dir, "Undetermined_S0_L00" + lane + "_I1_001.fastq.gz")
            unexpected_barcode_counts_i1 = re.sub(".fastq.gz", ".counts.txt", unaligned_i1)
            unaligned_i2 = ""
            outputs.append(unaligned_i1)
            if self.is_dual_index[lane]:
                unaligned_i2 = os.path.join(output_dir, "Undetermined_S0_L00" + lane + "_I2_001.fastq.gz")
                unexpected_barcode_counts_i2 = re.sub(".fastq.gz", ".counts.txt", unaligned_i2)
                outputs.append(unaligned_i2)
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
                        self.awk_read_1_processing_command(lane)
                    ),
                    bash.gzip(
                        None,
                        os.path.join(output_dir, "Undetermined_S0_L00" + lane + "_R1_001.fastq.gz")
                    ) if self.is_paired_end[lane] else None
                ],
                output_dependency=outputs,
                name= "fastq_convert.R1.unmatched." + self.run_id + "." + lane,
                samples=self.samples[lane]
            )
        )

        if self.is_paired_end[lane]:
            unmatched_R2_fastq = os.path.join(output_dir, "tmp", "unmatched_R2.fastq.gz")
            unaligned_i1 = os.path.join(output_dir, "Undetermined_S0_L00" + lane + "_I1_001.fastq.gz")
            unaligned_i2 = ""
            if unmatched_R2_fastq not in demuxfastqs_outputs:
                demuxfastqs_outputs.append(unmatched_R2_fastq)
            outputs = [
                os.path.join(output_dir, "Undetermined_S0_L00" + lane + "_R2_001.fastq.gz"),
                unaligned_i1
            ]
            unexpected_barcode_counts_i1 = re.sub(".fastq.gz", ".counts.txt", unaligned_i1)
#            demuxfastqs_outputs.append(unexpected_barcode_counts_i1)
            if self.is_dual_index[lane]:
                unaligned_i2 = os.path.join(output_dir, "Undetermined_S0_L00" + lane + "_I2_001.fastq.gz")
                outputs.append(unaligned_i2)
                unexpected_barcode_counts_i2 = re.sub(".fastq.gz", ".counts.txt", unaligned_i2)
#                demuxfastqs_outputs.append(unexpected_barcode_counts_i2)

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
                            self.awk_read_2_processing_command(lane)
                        )
                    ],
                    output_dependency=outputs,
                    name="fastq_convert.R2.unmatched." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
            )
        if unaligned_i1:
            postprocessing_jobs.append(
                run_processing_tools.fastq_unexpected_count(
                    unaligned_i1,
                    unexpected_barcode_counts_i1,
                    name="fastq_countbarcodes.I1.unmatched." + self.run_id + "." + lane,
                )
            )
        if unaligned_i2:
            postprocessing_jobs.append(
                run_processing_tools.fastq_unexpected_count(
                    unaligned_i2,
                    unexpected_barcode_counts_i2,
                    name="fastq_countbarcodes.I2.unmatched." + self.run_id + "." + lane,
                )
            )

        return demuxfastqs_outputs, postprocessing_jobs

    def awk_read_1_processing_command(self, lane, readset=None):
        """
        Returns a string serving as instructions for awk.
        This produces the command to convert the header of R1 fastq file from MGI to Illumina format
        For single-end libraies, the built command also extracts I1 (and I2 if exists) sequence from the R1 fastq,
        creating R1 (without barcode sequences) I1 (and I2) fastq files from R1 fastq file.
        """
        if self.is_paired_end[lane]:
            return """-v inst=\"{instrument}\" -v run=\"{run}\" 'match($0, /@({flowcell})L([0-9]{{1}})C([0-9]{{3}})R([0-9]{{3}})([0-9]+):?([ACTGN:-]+)?\/([0-9]{{1}})/, head_items) {{
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
        else:
            if self.is_dual_index[lane]:
                return """-v inst=\"{instrument}\" -v run=\"{run}\" -v read_len=\"{read_len}\" -v barcode1_len=\"{barcode1_len}\" -v barcode2_len=\"{barcode2_len}\" '{{
 header=$0
 getline seq
 getline sep
 getline qual
 match(header, /@({flowcell})L([0-9]{{1}})C([0-9]{{3}})R([0-9]{{3}})([0-9]+):?([ACTGN:-]+)?\/([0-9]{{1}})/, head_items)
 gsub("^0*", "", head_items[3]); gsub("^0*", "", head_items[4]); gsub("^0*", "", head_items[5])
 header="@" inst ":" run ":" head_items[1] ":" head_items[2] ":" head_items[5] ":" head_items[3] ":" head_items[4] " " head_items[7] ":N:0:" head_items[6]
 r1_seq=substr(seq,1,read_len)
 i1_seq=substr(seq,read_len+1,barcode1_len)
 i2_seq=substr(seq,read_len+barcode1_len+1,barcode2_len)
 r1_qual=substr(qual,1,read_len)
 i1_qual=substr(qual,read_len+1,barcode1_len)
 i2_qual=substr(qual,read_len+barcode1_len+1,barcode2_len)
 print header "\\n" r1_seq "\\n" sep "\\n" r1_qual | "gzip > {r1_out}"
 print header "\\n" i1_seq "\\n" sep "\\n" i1_qual | "gzip > {i1_out}"
 print header "\\n" i2_seq "\\n" sep "\\n" i2_qual | "gzip > {i2_out}"
}}'""".format(
                    instrument=self.instrument,
                    run=self.run_number,
                    read_len=self.read1cycles[lane],
                    barcode1_len=self.index1cycles[lane],
                    barcode2_len=self.index2cycles[lane],
                    flowcell=self.flowcell_id,
                    r1_out=readset.fastq1 if readset else os.path.join(self.output_dir, "Unaligned." + lane, "Undetermined_S0_L00" + lane + "_R1_001.fastq.gz"),
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
 r1_seq=substr(seq,1,read_len)
 i1_seq=substr(seq,read_len+1,barcode_len)
 r1_qual=substr(qual,1,read_len)
 i1_qual=substr(qual,read_len+1,barcode_len)
 print header "\\n" r1_seq "\\n" sep "\\n" r1_qual | "gzip > {r1_out}"
 print header "\\n" i1_seq "\\n" sep "\\n" i1_qual | "gzip > {i1_out}"
}}'""".format(
                    instrument=self.instrument,
                    run=self.run_number,
                    read_len=self.read1cycles[lane],
                    barcode_len=self.index1cycles[lane],
                    flowcell=self.flowcell_id,
                    r1_out=readset.fastq1 if readset else os.path.join(self.output_dir, "Unaligned." + lane, "Undetermined_S0_L00" + lane + "_R1_001.fastq.gz"),
                    i1_out=readset.index_fastq1 if readset else os.path.join(self.output_dir, "Unaligned." + lane, "Undetermined_S0_L00" + lane + "_I1_001.fastq.gz")
                )

    def awk_read_2_processing_command(self, lane, readset=None):
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
 match(header, /@({flowcell})L([0-9]{{1}})C([0-9]{{3}})R([0-9]{{3}})([0-9]+):?([ACTGN:-]+)?\/([0-9]{{1}})/, head_items)
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

    def throttle_jobs(self, jobs, key=""):
        """
        Group jobs of the same task (same name prefix) if they exceed the configured threshold number.
        """
        max_jobs_per_step = config.param('default', 'max_jobs_per_step', required=False, param_type="int")
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
                                name=f"{job_name}.{key}{('','.')[key!='']}{str(x + 1)}.{self.run_id}"
                            )
                        )
                reply.extend(merged_jobs)
            else:
                reply.extend(current_jobs)
        return reply

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
        if self.args.type == 'illumina':
            return parse_illumina_raw_readset_files(
                self.output_dir,
                self.run_dir,
                "PAIRED_END" if self.is_paired_end[lane] else "SINGLE_END",
                self.readset_file,
                lane,
                config.param('DEFAULT', 'genome_root', param_type="dirpath"),
                self.get_sequencer_minimum_read_length(),
                self.index1cycles[lane],
                self.index2cycles[lane],
                self.seqtype
            )
        else:
            return parse_mgi_raw_readset_files(
                self.readset_file,
                "PAIRED_END" if self.is_paired_end[lane] else "SINGLE_END",
                self.seqtype,
                self.run_id,
                lane,
                int(self.read1cycles[lane]),
                int(self.read2cycles[lane]),
                int(self.index1cycles[lane]),
                int(self.index2cycles[lane]),
                self.output_dir
            )

    def submit_jobs(self):
        super(RunProcessing, self).submit_jobs()

    @property
    def steps(self):
        # [ Illumina, mgig400, mgit7 ]
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
        utils.container_wrapper_argparse(argv)
    else:
        RunProcessing(protocol=['illumina', 'mgig400', 'mgit7'])

