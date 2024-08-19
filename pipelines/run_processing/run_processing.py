#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec
# Innovation Centre
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
from collections import Counter, OrderedDict

# Append genpipes directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs

from core.readset import parse_clarity_readset_file, parse_freezeman_readset_file
from bfx import bvatools
from bfx import picard
from bfx import fastp
from bfx import fastqc
from bfx import tools
from bfx import run_processing_tools
from bfx import bash_cmd as bash
from bfx import ngscheckmate
from bfx import bamixchecker
from bfx import gatk
from bfx.sequence_dictionary import split_by_size, parse_sequence_dictionary_file
from bfx import sambamba
from bfx import samtools

import utils

from pipelines import common

log = logging.getLogger(__name__)

class RunInfoRead(object):
    """
    Model of a read from the Illumina sequencer.
    Those attributes can be found in the RunInfo.xml file.
    """

    def __init__(self, number, nb_cycles, is_index, is_reverse_complement):
        self._number = number
        self._nb_cycles = nb_cycles
        self._is_index = is_index
        self._is_reverse_complement = is_reverse_complement

    @property
    def number(self):
        return self._number

    @property
    def nb_cycles(self):
        return self._nb_cycles

    @property
    def is_index(self):
        return self._is_index

    @property
    def is_reverse_complement(self):
        """
        New read feature from the release of RunInfo.xml v6.0, it is only
        leveraged by the NovaseqX.
        """
        return self._is_reverse_complement

class RunProcessing(common.MUGQICPipeline):
    """
    MGI Run Processing Pipeline
    ================================

    The standard Run Processing pipeline handles both Illumina and MGI
    sequencing technologies.  It uses the Illumina bcl2fastq software to
    convert and demultiplex Illumina base call files to fastq files.  In the
    case of MGI run processing, it uses fastq files produced by the MGI-G400
    sequencer, or MGI-T7 base call files, then does demultiplexing. Finally,
    the pipeline runs some QCs on the raw data, on the fastq and on the
    alignment.

    Sample Sheets
    -------------

    The pipeline uses one input sample sheet, a tsv file having the following
    columns:

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

        self.argparser.add_argument(
                "-t", "--type",
                help="""Sequencing technology : Illumina, MGI G400 or MGI T7
                (mandatory)""",
                choices=['illumina', 'mgig400', 'mgit7'],
                required=False)
        self.argparser.add_argument(
                "-r", "--readsets",
                help="Sample sheet for the MGI run to process (mandatory)",
                type=argparse.FileType('r'),
                required=False)
        self.argparser.add_argument(
                "-d", "--run",
                help="Run directory (mandatory)",
                required=False,
                dest="run_dir")
        self.argparser.add_argument(
                "--run-id",
                help="Run ID. Default is parsed from the run folder",
                required=False,
                dest="run_id")
        self.argparser.add_argument(
                "-f", "--flag",
                help="T7 flag files directory (mandatory for MGI T7 runs)",
                type=pathlib.Path, dest="raw_flag_dir",
                required=False)
        self.argparser.add_argument(
                "--splitbarcode-demux",
                help="""demultiplexing done while basecalling with MGI
                splitBarcode (only affect MGI G400 or T7 runs)""",
                action="store_true",
                required=False,
                dest="splitbarcode_demux")
        self.argparser.add_argument(
                "--lane",
                help="Lane number (to only process the given lane)",
                type=int,
                required=False,
                dest="lane_number")
        self.argparser.add_argument(
                "-x",
                help="""First index base to use for demultiplexing (inclusive).
                The index from the sample sheet will be adjusted according to
                that value.""",
                type=int,
                required=False,
                dest="first_index")
        self.argparser.add_argument(
                "-y",
                help="Last index base to use for demultiplexing (inclusive)",
                type=int,
                required=False,
                dest="last_index")
        self.argparser.add_argument(
                "-m",
                help="""Number of index mistmaches allowed for demultiplexing
                (default 1). Barcode collisions are always checked.""",
                type=int,
                required=False,
                dest="number_of_mismatches")
        self.argparser.add_argument(
                "--allow-barcode-collision",
                help="""Allow barcode collision by not comparing barcode
                sequences to each other (usually decreases the demultiplexing
                efficiency).""",
                action="store_true",
                required=False,
                dest="allow_barcode_collision")

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
            _raise(SanitycheckError(f"Unsupported protocol {typearg}"))
        if flag and not typearg == 'mgit7':
            log.info("Ignoring -f/--flag option because useless without '-t/--type mgit7'...")
        if typearg == 'illumina' and splitbarcode:
            log.info("Ignoring --splitbarcode-demux option because useless with '-t/--type illunina'...")

        super(RunProcessing, self).__init__(protocol)

    @property
    def output_dirs(self):
        if not hasattr(self, "_output_dirs"):
            self._output_dirs = {}
            for lane in self.lanes:
                self._output_dirs[lane] = {
                    f"Unaligned.{lane}_directory": os.path.relpath(os.path.join(self.output_dir, f"Unaligned.{lane}"), self.output_dir),
                    f"Unaligned.{lane}_directory_tmp": os.path.relpath(os.path.join(self.output_dir, f"Unaligned.{lane}_tmp"), self.output_dir),
                    f"Unaligned.{lane}.noindex_directory": os.path.relpath(os.path.join(self.output_dir, f"Unaligned.{lane}.noindex"), self.output_dir),
                    f"Aligned.{lane}_directory": os.path.relpath(os.path.join(self.output_dir, f"Aligned.{lane}"), self.output_dir),
                    "index_directory": os.path.relpath(os.path.join(self.output_dir, 'index'), self.output_dir),
                    "report_directory": os.path.relpath(os.path.join(self.output_dir, 'report'), self.output_dir)
                }
        return self._output_dirs

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
        For Illumina runs :
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
                        _raise(SanitycheckError(f"Error: Run ID could not be parsed from the RUN folder : {self.run_dir}"))
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
                _raise(SanitycheckError(f"Error: Multiple flowcells found in this run. Please check the readset file. {msg}"))
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
                if os.path.exists(os.path.join(self.output_dirs[lane][f"Unaligned.{lane}_directory"], "raw_fastq", "BioInfo.csv")):
                    self._bioinfo_files[lane] = os.path.join(self.output_dirs[lane][f"Unaligned.{lane}_directory"], "raw_fastq", "BioInfo.csv")
                else:
                    bioinfo_file = os.path.join(self.run_dir, f"L0{lane}", "BioInfo.csv")
                    if lane in self._bioinfo_files and self._bioinfo_files[lane] != bioinfo_file:
                        _raise(SanitycheckError(f"More than one Bioinfo.csv found for lane '{lane}' : {self._bioinfo_files[lane]}, {bioinfo_file}"))
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
            self._json_flag_files = {}
            for lane in self.lanes:
                json_flag_file = os.path.join(self.output_dir, f"flag.{lane}.json")
                for filename in os.listdir(self.raw_flag_dir):
                    if re.match(f"{self.run_id}_{lane}_.+json", filename):
                        if not os.path.exists(json_flag_file):
                            if not os.path.exists(os.path.dirname(json_flag_file)):
                                os.makedirs(os.path.dirname(json_flag_file))
                            shutil.copy(os.path.join(self.raw_flag_dir, filename), json_flag_file)
                        self._json_flag_files[lane] = json_flag_file
                        log.info(f"JSON FLAG file for lane {lane} : {json_flag_file}")
                        break
                else:
                    _raise(SanitycheckError(f"Could not find any proper JSON flag file in {self.raw_flag_dir} for RUN {self.run_id}"))
        return self._json_flag_files

    @property
    def barcode_files(self):
        if not hasattr(self, "_barcode_files"):
            self._barcode_files = {}
            for lane in self.lanes:
                barcode_file = os.path.join(self.output_dir, f"barcodes.{lane}.txt")
                if not os.path.exists(barcode_file):
                    if not os.path.exists(os.path.dirname(barcode_file)):
                        os.makedirs(os.path.dirname(barcode_file))
                self._barcode_files[lane] = barcode_file
                log.info(f"BARCODE file for lane {lane} : {barcode_file}")
        return self._barcode_files

    @property
    def json_flag_hash(self):
        if not hasattr(self, "_json_flag_hash"):
            self._json_flag_hash = {}
            for lane in self.lanes:
                with open(self.json_flag_files[lane], "r") as jff:
                    json_flag_content = json.load(jff)
                # Turns the InfoVector into a dictionnary and concat it back to
                # the flagfile json dict
                keys = [item[0] for item in json_flag_content['experimentInfoVec']]
                vals = [item[1] for item in json_flag_content['experimentInfoVec']]
                del json_flag_content['experimentInfoVec']
                self._json_flag_hash[lane] = dict(zip(keys, vals)) | json_flag_content
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
                if is_json(self.readset_file):
                    self._lanes = [str(lane) for lane in list(set([sample['lane'] for sample in json.load(open(self.readset_file, 'r'))['samples']]))]
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
        return self._index1cycles

    @property
    def index2cycles(self):
        if not hasattr(self, "_index2cycles"):
            self._index2cycles = {}
            for lane in self.lanes:
                self._index2cycles[lane] = self.get_index2cycles(lane)
        return self._index2cycles

    @property
    def index1reverse(self):
        if not hasattr(self, "_index1reverse"):
            self._index1reverse = {}
            for lane in self.lanes:
                self._index1reverse[lane] = self.get_indexreverse(lane)[0]
        return self._index1reverse

    @property
    def index2reverse(self):
        if not hasattr(self, "_index2reverse"):
            self._index2reverse = {}
            for lane in self.lanes:
                self._index2reverse[lane] = self.get_indexreverse(lane)[1]
        return self._index2reverse

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
    def read_infos(self):
        if not hasattr(self, "_read_infos"):
            self._read_infos = self.parse_run_info_file()
        return self._read_infos

    @property
    def sbs_consumable_version(self):
        # For Illumina runs only
        if not hasattr(self, "_sbs_consumable_version"):
            self._sbs_consumable_version = self.get_sbs_consumable_version()
        return self._sbs_consumable_version

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
                    _raise(SanitycheckError(f"One run (\"{self.run_id}\") cannot be shared by different instruments !! (\"{lane_instrument}\" vs. \"{self._instrument}\")"))
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
                    _raise(SanitycheckError(f"One run (\"{self.run_id}\") cannot be placed on multiple flowcell positions !! (\"{lane_flowcell_position}\" vs. \"{self._flowcell_position}\")"))
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
                    _raise(SanitycheckError(f"One run (\"{self.run_id}\") cannot be defined by more than one run counter !! (\"{lane_run_number}\" vs. \"{run_number}\")"))
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
                    _raise(SanitycheckError(f"Sequencer Run ID conflct (\"{lane_sequencer_run_id}\" vs. \"{sequencer_run_id}\")"))
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
        Get year of the run from sample sheet
        """
        if not hasattr(self, "_year"):
            if is_json(self.readset_file):
                dates = [json.load(open(self.readset_file, 'r'))['run_start_date']]
            else:
                dates = set([date for date in list(set([line['Start Date'] for line in csv.DictReader(open(self.readset_file, 'r'), delimiter='\t', quotechar='"')]))])
            if len(list(dates)) > 1:
                _raise(SanitycheckError(f"More than one date were found in the sample sheet for the run \"{self.run_id}\""))
            else:
                self._year = list(dates)[0].split("-")[0]
        return self._year

    @property
    def date(self):
        """
        Get whole date of the run from sample sheet
        """
        if not hasattr(self, "_date"):
            if is_json(self.readset_file):
                dates = [json.load(open(self.readset_file, 'r'))['run_start_date']]
            else:
                dates = set([date for date in list(set([line['Start Date'] for line in csv.DictReader(open(self.readset_file, 'r'), delimiter='\t', quotechar='"')]))])
            if len(list(dates)) > 1:
                _raise(SanitycheckError(f"More than one date were found in the sample sheet for the run \"{self.run_id}\""))
            else:
                date = list(dates)[0].split("-")
                self._date = date[0][-2:] + date[1] + date[2]
        return self._date

    @property
    def report_hash(self):
        if not hasattr(self, "_report_hash"):
            self._report_hash = {}


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
                    os.path.basename(self.run_dir.rstrip('/')) + "-" + self.seqtype
                    )
            for lane in self.lanes:
                self._report_hash[lane] = {
                    "version" : "3.0",
                    "run" : self.run_id,
                    "run_obj_id": self.readsets[lane][0].run_obj_id if self.readsets[lane][0].run_obj_id else None,
                    "instrument" : self.instrument,
                    "flowcell" : self.flowcell_id,
                    "lane" : lane,
                    "seqtype" : self.seqtype,
                    "sequencing_method" : "PAIRED_END" if self.is_paired_end[lane] else "SINGLE_END",
                    "steps" : [],
                    "readsets" : dict(
                        [
                            (
                                readset.name,
                                {
                                    "sample_name": readset.sample.name,
                                    "barcodes": readset.indexes,
                                    "species": readset.species,
                                    "reported_sex": readset.gender,
                                    "pool_fraction": readset.pool_fraction,
                                    "library_type": readset.protocol,
                                    "fastq_1": {
                                            "path": readset.fastq1,
                                            "size": None,
                                            "final_path": os.path.join(full_destination_folder, readset.fastq1)
                                            },
                                    "fastq_2": {
                                            "path": readset.fastq2 if self.is_paired_end[lane] else None,
                                            "size": None,
                                            "final_path": os.path.join(full_destination_folder, readset.fastq2) if self.is_paired_end[lane] else None
                                            },
                                    "bam": {
                                            "path": readset.bam + ".bam" if readset.bam else None,
                                            "size": None,
                                            "final_path": os.path.join(full_destination_folder, readset.bam + ".bam") if readset.bam else None
                                            },
                                    "bai": {
                                            "path": readset.bam + ".bai" if readset.bam else None,
                                            "size": None,
                                            "final_path": os.path.join(full_destination_folder, readset.bam + ".bai") if readset.bam else None
                                            },
                                    "derived_sample_obj_id": readset.library,
                                    "project_obj_id": readset.project_id,
                                    "project_name": readset.project_name,
                                    "external_project_id": readset.external_project_id if is_json(self.readset_file) else None
                                }
                            ) for readset in self.readsets[lane]
                        ]
                    )
                }
        return self._report_hash

    @property
    def run_validation_report_json(self):
        if not hasattr(self, "_run_validation_report_json"):
            self._run_validation_report_json = {}
            for lane in self.lanes:
                self._run_validation_report_json[lane] = os.path.join(self.output_dirs[lane]["report_directory"], f"{self.run_id}.{lane}.run_validation_report.json")
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

            unaligned_dir = self.output_dirs[lane][f"Unaligned.{lane}_directory"]
            basecall_dir = os.path.join(unaligned_dir, "basecall")

            lane_config_file = os.path.join(unaligned_dir, f"{self.run_id}.{lane}.settings.config")

            # If demultiplexing is perform while basecalling
            if self.args.splitbarcode_demux:
                # Add the barcodes in the JSON flag file
                self.edit_mgi_t7_flag_file(lane)
                self.create_barcode_file(lane)

                basecall_outputs, postprocessing_jobs = self.generate_basecall_outputs(lane)
                basecall_outputs.extend(
                    [
                        os.path.join(basecall_dir, self.run_id, f"L0{lane}"),
                        os.path.join(basecall_dir, self.run_id, f"L0{lane}", "SequenceStat.txt"),
                        os.path.join(basecall_dir, self.run_id, f"L0{lane}", "BarcodeStat.txt")
                    ]
                )

                if self.run_id == self.flowcell_id:
                    symlink_command = None
                else:
                    symlink_command = bash.ln(
                            os.path.relpath(os.path.join(basecall_dir, self.flowcell_id), basecall_dir),
                            os.path.join(basecall_dir, self.run_id),
                            input = os.path.join(basecall_dir, self.flowcell_id)
                            )

                lane_jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(basecall_dir),
                            run_processing_tools.mgi_splitbarcode(
                                input,
                                self.run_dir,
                                self.flowcell_id,
                                basecall_outputs,
                                basecall_dir,
                                self.json_flag_hash[lane],
                                self.barcode_files[lane],
                                self.number_of_mismatches
                            ),
                            symlink_command
                        ],
                        name=f"basecall.{self.run_id}.{lane}",
                        samples=self.samples[lane],
                        input_dependency=[input],
                        report_files=[os.path.join(basecall_dir, self.run_id, f"L0{lane}", "SequenceStat.txt")]
                    )
                )
                for readset in self.readsets[lane]:
                    readset.report_files['basecall'] = [os.path.join(basecall_dir, self.run_id, f"L0{lane}", "SequenceStat.txt")]

                if postprocessing_jobs:
                    jobs_to_throttle.extend(postprocessing_jobs)

                # parse splitBarcode metrics to be compatible with multiqc
                parse_job = concat_jobs(
                        [
                            run_processing_tools.parse_splitBarcode_metrics(
                                os.path.join(basecall_dir, self.run_id, f"L0{lane}", "BarcodeStat.txt"),
                                os.path.join(self.output_dir, "samplesheet." + lane + ".csv"),
                                os.path.join(basecall_dir, self.run_id, f"L0{lane}", "BarcodeStat_L0" + lane + "_multiqc.txt")
                                )
                        ],
                        name=f"parse_splitBarcode_metrics.{self.run_id}.{lane}",
                        input_dependency=[os.path.join(basecall_dir, self.run_id, f"L0{lane}", "BarcodeStat.txt")],
                        samples=self.samples[lane],
                        report_files=[os.path.join(basecall_dir, self.run_id, f"L0{lane}", "BarcodeStat_L01" + lane + "_multiqc.txt")]
                    )

                lane_jobs.append(parse_job)

            else:
                basecall_outputs = [
                    os.path.join(basecall_dir, self.run_id, f"L0{lane}"),
                    os.path.join(basecall_dir, self.run_id, f"L0{lane}", f"{self.raw_fastq_prefix}_L0{lane}_read_1.fq.gz"),
                    os.path.join(basecall_dir, self.run_id, f"L0{lane}", f"{self.raw_fastq_prefix}_L0{lane}_read_2.fq.gz"),
                    os.path.join(basecall_dir, self.run_id, f"L0{lane}", f"{self.raw_fastq_prefix}_L0{lane}.summaryReport.html"),
                    os.path.join(basecall_dir, self.run_id, f"L0{lane}", f"{self.raw_fastq_prefix}_L0{lane}.heatmapReport.html"),
                    os.path.join(basecall_dir, self.run_id, f"L0{lane}", "summaryTable.csv")
                ]

                raw_fastq_dir = os.path.join(unaligned_dir, "raw_fastq")
                raw_fastq_outputs = [
                    os.path.join(raw_fastq_dir, f"{self.raw_fastq_prefix}_L0{lane}_read_1.fq.gz"),
                    os.path.join(raw_fastq_dir, f"{self.raw_fastq_prefix}_L0{lane}_read_2.fq.gz"),
                    os.path.join(raw_fastq_dir, f"{self.raw_fastq_prefix}_L0{lane}.summaryReport.html"),
                    os.path.join(raw_fastq_dir, f"{self.raw_fastq_prefix}_L0{lane}.heatmapReport.html"),
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
                            os.path.relpath(os.path.join(basecall_dir, self.run_id, f"L0{lane}"), os.path.dirname(raw_fastq_dir)),
                            raw_fastq_dir,
                            input=os.path.join(basecall_dir, self.run_id, f"L0{lane}")
                        )
                    ],
                    name=f"basecall.{self.run_id}.{lane}",
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
        Generate a file with all the indexes found in the index-reads of the
        run.

        The input barcode file is a two columns tsv file. Each line has a
        `barcode_sequence` and the corresponding `barcode_name`. This file can
        be generated by a LIMS.

        The output is a tsv file named `RUNFOLDER_LANENUMBER.metrics` that will
        be saved in the output directory. This file has four columns, the
        barcode/index sequence, the index name, the number of reads and the
        number of reads that have passed the filter.
        """
        jobs = []

        # Instantiate the readset to have the readset file a dependency
        # and self.mask available for NovaseqX
        self.readsets

        # Builds the mask from the ReadInfo.xml
        cycles_mask = ""
        for read in self.read_infos:
            if read.is_index:
                cycles_mask += str(read.nb_cycles) + "B"
                break
            else:
                cycles_mask += str(read.nb_cycles) + 'T'
        if "B" not in cycles_mask:
            log.info(" ".join(["No indexes cycles set on the sequencer,",
                               "*NOT* Generating index counts"]))
        else:
            for lane in self.lanes:
                if int(self.index1cycles[lane]) + int(self.index2cycles[lane]) == 0 and len(self.readsets[lane]) > 1:
                    err_msg = "LANE SETTING ERROR :\n"
                    err_msg += f"Unable to demultiplex {str(len(self.readsets[lane]))} samples : No barcode in fastq files...\n"
                    err_msg += f"(in {self.run_dir})"
                    _raise(SanitycheckError(err_msg))

                lane_jobs = []

                # If masks are defined per lanes in the RunProcessing object
                # use them .
                if self.mask[lane]:
                    mask = ""
                    for component in self.mask[lane].split(","):
                        if component[0] == "Y":
                            mask += re.search("\d+",component)[0] + "T"
                        elif component[0] == "I":
                            mask += re.search("\d+",component)[0] + "B"
                            break
                else:
                    mask = cycles_mask

                input = os.path.join(self.run_dir, "RunInfo.xml")
                output = os.path.join(self.output_dirs[lane]["index_directory"],  f"{self.run_id}_{lane}.metrics")
                basecalls_dir = os.path.join(self.output_dirs[lane]["index_directory"], "BaseCalls")

                barcode_file = config.param('index', 'barcode_file', param_type='filepath', required=False)
                if not (barcode_file and os.path.isfile(barcode_file)):
                    resources_dir = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'resources')
                    if (self.seqtype in ["hiseqx", "hiseq4000", "iSeq"] \
                            or (self.seqtype in ["novaseq"] \
                            and self.sbs_consumable_version == 3)):
                        barcode_file = os.path.join(resources_dir, 'barcodes_by_sequence.i5rev.txt')
                    else:
                        barcode_file = os.path.join(resources_dir, 'barcodes_by_sequence.i5fwd.txt')

                # CountIlluminaBarcode
                lane_jobs.append(
                    concat_jobs(
                    [
                            bash.mkdir(basecalls_dir),
                            bash.rm(
                                os.path.join(basecalls_dir, f"L00{lane}"),
                                force=True
                            ),
                            bash.ln(
                                os.path.join(self.run_dir, "Data", "Intensities", "BaseCalls", f"L00{lane}"),
                                os.path.join(basecalls_dir, f"L00{lane}")
                            ),
                            bash.ln(
                                os.path.join(self.run_dir, "Data", "Intensities", "s.locs"),
                                os.path.join(self.output_dirs[lane]["index_directory"], "s.locs")
                            ),
                            run_processing_tools.index(
                                input,
                                barcode_file,
                                basecalls_dir,
                                self.number_of_mismatches,
                                lane,
                                mask,
                                output
                            )
                        ],
                        name=f"index.{self.run_id}.{lane}",
                        samples=self.samples[lane],
                        input_dependency=[self.readset_file,
                                          os.path.join(self.run_dir,
                                                       "RunInfo.xml")],
                        report_files=[output]
                    )
                )

                # for readset in self.readsets[lane]:
                #     readset.report_files['index'] = [output]

                self.add_to_report_hash("index", lane,  lane_jobs)
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
            jobs_to_throttle = []

            input = self.readset_file
            fastq_outputs, postprocessing_jobs, cleanjob_deps = self.generate_bcl2fastq_outputs(lane)
            output_dir = self.output_dirs[lane][f"Unaligned.{lane}_directory"]
            tmp_output_dir = f"{output_dir}_tmp"
            casava_sample_sheet = os.path.join(self.output_dir, f"casavasheet.{lane}.indexed.csv")

            fastq_outputs.append(os.path.join(tmp_output_dir, "Stats/Stats.json"))

            demultiplexing_done_file = os.path.join(self.output_dir, f"{self.run_id}.{lane}.demultiplexingDone")

            if os.path.exists(demultiplexing_done_file) and not self.force_jobs:
                log.info(f"Demultiplexing done already... Skipping fastq step for lane {lane}...")

            else:

                if self.umi:
                    output_dir_noindex = self.output_dirs[lane][f"Unaligned.{lane}.noindex_directory"]
                    casava_sample_sheet_noindex = os.path.join(self.output_dir, f"casavasheet.{lane}.noindex.csv")

                    lane_jobs.append(
                        concat_jobs(
                            [
                                run_processing_tools.bcl2fastq(
                                    input,
                                    fastq_outputs,
                                    tmp_output_dir,
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
                                )
                            ],
                            name=f"fastq_illumina.{self.run_id}.{lane}",
                            samples=self.samples[lane]
                        )
                    )

                else:
                    bcl2fastq_job = run_processing_tools.bcl2fastq(
                        input,
                        fastq_outputs,
                        tmp_output_dir,
                        casava_sample_sheet,
                        self.run_dir,
                        lane,
                        self.bcl2fastq_extra_option,
                        demultiplex=True,
                        mismatches=self.number_of_mismatches,
                        mask=self.mask[lane],
                        ini_section='fastq_illumina'
                    )
                    bcl2fastq_job.name = f"fastq_illumina.{self.run_id}.{lane}"
                    bcl2fastq_job.samples = self.samples[lane]
                    lane_jobs.append(bcl2fastq_job)

                postprocessing_jobs.append(
                    concat_jobs(
                        [
                            bash.mv(os.path.join(tmp_output_dir, "Stats"), output_dir, force=True),
                            bash.mv(os.path.join(tmp_output_dir, "Reports"), output_dir, force=True),
                            bash.mv(os.path.join(tmp_output_dir, f"Undetermined_S0_L00{lane}_R1_001.fastq.gz"), output_dir, force=True),
                            bash.mv(os.path.join(tmp_output_dir, f"Undetermined_S0_L00{lane}_I1_001.fastq.gz"), output_dir, force=True),
                            bash.mv(os.path.join(tmp_output_dir, f"Undetermined_S0_L00{lane}_R2_001.fastq.gz"), output_dir, force=True) if self.is_paired_end[lane] else None,
                            bash.mv(os.path.join(tmp_output_dir, f"Undetermined_S0_L00{lane}_I2_001.fastq.gz"), output_dir, force=True) if self.is_dual_index[lane] else None
                        ],
                        report_files=[os.path.join(output_dir, "Stats/Stats.json")],
                        samples=self.samples[lane],
                        name=f"fastq_move.undet_reports_stats.{self.run_id}.{lane}",
                        output_dependency=[os.path.join(output_dir, "Stats/Stats.json")]
                    )
                )

                for readset in self.readsets[lane]:
                    readset.report_files['fastq'] = [os.path.join(output_dir, "Stats/Stats.json")]

                # Last post-processing job : cleaning of the tmp folder
                postprocessing_jobs.append(
                    concat_jobs(
                        [
                            bash.rm(
                                tmp_output_dir,
                                recursive=True,
                                force=True
                            ),
                            bash.touch(demultiplexing_done_file)
                        ],
                        input_dependency=cleanjob_deps,
                        name=f"fastq_clean_tmp.{self.run_id}.{lane}",
                        samples=self.samples[lane]
                    )
                )

                jobs_to_throttle.extend(postprocessing_jobs)
                lane_jobs.extend(self.throttle_jobs(jobs_to_throttle, lane))

                # don't depend on notification commands
                self.add_copy_job_inputs(lane_jobs, lane)

                notification_command_start = config.param('fastq_notification_start', 'notification_command', required=False)
                if notification_command_start:
                    notification_command_start = notification_command_start.format(
                        output_dir=self.output_dir,
    #                    number_of_mismatches=self.number_of_mismatches,
                        lane_number=lane,
    #                    mask=self.mask[lane],
    #                    technology=config.param('fastq', 'technology'),
    #                    run_id=self.run_id
                    )
                    # Use the same inputs and output of fastq job to send a notification each time the fastq job run
                    job = Job(
                        [input],
                        [f"notificationFastqStart.{lane}.out"],
                        command=notification_command_start,
                        name=f"fastq_notification_start.{self.run_id}.{lane}",
                        samples=self.samples[lane]
                    )
    #                lane_jobs.append(job)

                notification_command_end = config.param('fastq_notification_end', 'notification_command', required=False)
                if notification_command_end:
                    notification_command_end = notification_command_end.format(
                        output_dir=self.output_dir,
                        lane_number=lane,
    #                    technology=config.param('fastq', 'technology'),
    #                    run_id=self.run_id
                    )
                    job = Job(
                        fastq_outputs,
                        [f"notificationFastqEnd.{lane}.out"],
                        command=notification_command_end,
                        name=f"fastq_notification_end.{self.run_id}.{lane}",
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
            _raise(SanitycheckError(f"Could not determine which section to use for fastq step from given protocol {self.protocol}"))
        jobs = []

        if self.args.splitbarcode_demux:
            log.info("Demultiplexing done during the basecalling... Skipping fastq step...")

        else:
            log.info("Start demultiplexing the FASTQ files...")

            for lane in self.lanes:

                demultiplexing_done_file = os.path.join(self.output_dir, f"{self.run_id}.{lane}.demultiplexingDone")

                if int(self.index1cycles[lane]) + int(self.index2cycles[lane]) == 0 and len(self.readsets[lane]) > 1:
                    err_msg = "LANE SETTING ERROR :\n"
                    err_msg += f"Unable to demultiplex {str(len(self.readsets[lane]))} samples : No barcode in fastq files...\n"
                    err_msg += f"(in {self.run_dir})"
                    _raise(SanitycheckError(err_msg))

                elif os.path.exists(demultiplexing_done_file) and not self.force_jobs:
                    log.info(f"Demultiplexing done already... Skipping fastq step for lane {lane}...")

                else:

                    lane_jobs = []
                    jobs_to_throttle = []

                    unaligned_dir = self.output_dirs[lane][f"Unaligned.{lane}_directory"]
                    raw_fastq_dir = os.path.join(unaligned_dir, "raw_fastq")
                    raw_name_prefix = f"{self.raw_fastq_prefix}_L0{lane}"

                    input1 = os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read_1.fq.gz")
                    input2 = os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read_2.fq.gz")

                    # Prepare the link job here, even though it will only be used in MGI-G400 cases
                    link_raw_fastq_job = concat_jobs(
                        [
                            bash.mkdir(raw_fastq_dir),
                            bash.ln(
                                os.path.join(self.run_dir, f"L0{lane}", "*"),
                                raw_fastq_dir,
                            )
                        ]
                    )
                    link_raw_fastq_job.output_files = [
                        os.path.join(raw_fastq_dir, f"{raw_name_prefix}.summaryReport.html"),
                        os.path.join(raw_fastq_dir, f"{raw_name_prefix}.heatmapReport.html"),
                        os.path.join(raw_fastq_dir, "summaryTable.csv")
                    ]
                    if self.readsets[lane][0].run_type == "PAIRED_END":
                        link_raw_fastq_job.output_files.append(input1)
                        link_raw_fastq_job.output_files.append(input2)
                        link_raw_fastq_job.output_files.append(os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read_1.fq.fqStat.txt"))
                        link_raw_fastq_job.output_files.append(os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read_2.fq.fqStat.txt"))
                    else:
                        link_raw_fastq_job.output_files.append(os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read.fq.gz"))
                        link_raw_fastq_job.output_files.append(os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read.fq.fqStat.txt"))
                    link_raw_fastq_job.output_files.append(os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read.report.html"))
                    link_raw_fastq_job.input_files = [re.sub(raw_fastq_dir, os.path.join(self.run_dir, f"L0{lane}"), dep) for dep in link_raw_fastq_job.output_files]

                    if (len(self.readsets[lane]) == 1) and (int(self.index1cycles[lane]) + int(self.index2cycles[lane]) == 0):
                        log.info(f"No barcode cycles in the lane... Fastq files will just be copied and renamed, for lane {lane}...")

                        readset = self.readsets[lane][0]
                        if readset.run_type == "PAIRED_END":
                            demultiplex_job = concat_jobs(
                                [
                                    link_raw_fastq_job if (ini_section == 'fastq_g400') else None,
                                    bash.mkdir(os.path.dirname(readset.fastq1)),
                                    bash.cp(
                                        os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read_1.fq.gz"),
                                        readset.fastq1,
                                    ),
                                    bash.cp(
                                        os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read_1.fq.fqStat.txt"),
                                        re.sub("gz", "fqStat.txt", readset.fastq1),
                                    ),
                                    bash.cp(
                                        os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read.report.html"),
                                        re.sub("_R1_001.fastq.gz", "_read.report.html", readset.fastq1),
                                    ),
                                    bash.cp(
                                        os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read_2.fq.gz"),
                                        readset.fastq2,
                                    ),
                                    bash.cp(
                                        os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read_2.fq.fqStat.txt"),
                                        re.sub("gz", "fqStat.txt", readset.fastq2),
                                    )
                                ]
                            )
                        else:
                            demultiplex_job = concat_jobs(
                                [
                                    link_raw_fastq_job if (ini_section == 'fastq_g400') else None,
                                    bash.mkdir(os.path.dirname(readset.fastq1)),
                                    bash.cp(
                                        os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read.fq.gz"),
                                        readset.fastq1,
                                    ),
                                    bash.cp(
                                        os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read.fq.fqStat.txt"),
                                        re.sub("", "fqStat.txt", readset.fastq1),
                                    ),
                                    bash.cp(
                                        os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read.report.html"),
                                        re.sub("_R1_001.fastq.gz", "_read.report.html", readset.fastq1),
                                    )
                                ]
                            )
                        demultiplex_job.name = f"{ini_section}.demultiplex.{self.run_id}.{lane}"
                        demultiplex_job.samples = self.samples[lane]
                        lane_jobs.append(demultiplex_job)

                    else:
                        # Demultiplexing
                        demuxfastqs_outputs, postprocessing_jobs, cleanjob_deps = self.generate_demuxfastqs_outputs(lane)

                        tmp_output_dir = os.path.dirname(demuxfastqs_outputs[0])
                        tmp_metrics_file = os.path.join(tmp_output_dir, self.run_id + "." + lane + ".DemuxFastqs.metrics.txt")
                        metrics_file = os.path.join(self.output_dirs[lane][f"Unaligned.{lane}_directory"], self.run_id + "." + lane + ".DemuxFastqs.metrics.txt")
                        demuxfastqs_outputs.append(metrics_file)

                        if self.readsets[lane][0].run_type == "PAIRED_END":
                            raw_name_prefix = f"{self.raw_fastq_prefix}_L0{lane}"
                            input1 = os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read_1.fq.gz")
                            input2 = os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read_2.fq.gz")
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
                                ini_section=ini_section
                            )
                        else:
                            input1 = os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read.fq.gz")
                            if ini_section == 'fastq_g400':
                                input1 = os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read.fq.gz")
                                demultiplex_job = run_processing_tools.demux_fastqs_single_end(
                                    os.path.join(self.output_dir, "samplesheet." + lane + ".csv"),
                                    self.number_of_mismatches,
                                    self.mask[lane],
                                    demuxfastqs_outputs,
                                    tmp_metrics_file,
                                    input1,
                                    ini_section=ini_section
                                )
                            elif ini_section == 'fastq_t7':
                                demultiplex_job = run_processing_tools.demux_fastqs_by_chunk(
                                    os.path.join(self.output_dir, "samplesheet." + lane + ".csv"),
                                    self.number_of_mismatches,
                                    self.mask[lane],
                                    demuxfastqs_outputs,
                                    tmp_metrics_file,
                                    input1,
                                    ini_section=ini_section
                                )
                        lane_jobs.append(
                            concat_jobs(
                                [
                                    link_raw_fastq_job if (ini_section == 'fastq_g400') else None,
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
                                name=f"{ini_section}.demultiplex.{self.run_id}.{lane}",
                                samples=self.samples[lane],
                                input_dependency=link_raw_fastq_job.input_files if (ini_section == 'fastq_g400') else demultiplex_job.input_files,
                                output_dependency=demuxfastqs_outputs + [tmp_output_dir] + link_raw_fastq_job.output_files if (ini_section == 'fastq_g400') else demuxfastqs_outputs + [tmp_output_dir],
                                report_files=[metrics_file]
                            )
                        )
                        for readset in self.readsets[lane]:
                            readset.report_files['fastq'] = [metrics_file]

                        # Last post-processing job : cleaning of the tmp folder
                        postprocessing_jobs.append(
                            concat_jobs(
                                [
                                    bash.rm(
                                        tmp_output_dir,
                                        recursive=True,
                                        force=True
                                    ),
                                    bash.touch(demultiplexing_done_file)
                                ],
                                input_dependency=cleanjob_deps,
                                name=f"fastq_clean_tmp.{self.run_id}.{lane}",
                                samples=self.samples[lane]
                            )
                        )

                        jobs_to_throttle.extend(postprocessing_jobs)
                        if self.args.type == 'mgit7':
                            lane_jobs.extend(jobs_to_throttle)
                        else:
                            lane_jobs.extend(self.throttle_jobs(jobs_to_throttle, lane))

                    self.add_to_report_hash("fastq", lane, lane_jobs)
                    self.add_copy_job_inputs(lane_jobs, lane)
                    jobs.extend(lane_jobs)

        return jobs

    ## DEPRECATED
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
                        name=f"qc_graphs.{readset.name}.qc.{self.run_id}.{lane}",
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
                                        readset.index_fastq2 if readset.index_type == "DUALINDEX" else None,
                                        "FASTQ",
                                        region_name,
                                        output_dir + "_index"
                                    )
                                ],
                                name=f"qc_graphs.{readset.name}.qc_index.{self.run_id}.{lane}",
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

    ## DEPRECATED
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
                            readset.report_files['fastqc'] = [outputs[2]]
                    elif input1 == readset.fastq2:
                        job_suffix = "R2."
                        input2 = None
                        unzip=True
                        if not len(self.readsets[lane]) == 1:
                            readset.report_files['fastqc'].append(outputs[2])
                    elif input1 == readset.index_fastq1:
                        job_suffix = "Barcodes."
                        input2 = readset.index_fastq2 if readset.index_type == "DUALINDEX" else None
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
                            name=f"fastqc.{readset.name}_{job_suffix}.{self.run_id}.{lane}",
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
            if 'mgi' in self.args.type and not self.args.splitbarcode_demux:
                unaligned_dir = self.output_dirs[lane][f"Unaligned.{lane}_directory"]
                raw_fastq_dir = os.path.join(unaligned_dir, "raw_fastq")
                raw_name_prefix = f"{self.raw_fastq_prefix}_L0{lane}"
                raw_reads_input1 = os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read_1.fq.gz")
                raw_reads_input2 = os.path.join(raw_fastq_dir, f"{raw_name_prefix}_read_2.fq.gz")
                raw_reads_output_json_path = os.path.join(raw_fastq_dir, "fastp", f"{raw_name_prefix}.raw_reads.fastp.json")
                raw_reads_output_html_path = os.path.join(raw_fastq_dir, "fastp", f"{raw_name_prefix}.raw_reads.fastp.html")
                lane_jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                os.path.dirname(raw_reads_output_json_path),
                                remove=True
                            ),
                            fastp.basic_qc(
                                raw_reads_input1,
                                raw_reads_input2,
                                raw_reads_output_json_path,
                                raw_reads_output_html_path
                            )
                        ],
                        name=f"fastp.raw_reads.{self.run_id}.{lane}",
                        samples=self.samples[lane],
                        output_dependency=[raw_reads_output_json_path, raw_reads_output_html_path],
                        report_files=[raw_reads_output_json_path]
                    )
                )

            for readset in self.readsets[lane]:
                output_json_path = os.path.join(os.path.dirname(readset.fastq1), "fastp", readset.name + ".fastp.json")
                output_html_path = os.path.join(os.path.dirname(readset.fastq1), "fastp", readset.name + ".fastp.html")
                lane_jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                os.path.dirname(output_json_path),
                                remove=True
                            ),
                            fastp.basic_qc(
                                readset.fastq1,
                                readset.fastq2,
                                output_json_path,
                                output_html_path
                            )
                        ],
                        name=f"fastp.{readset.name}.{self.run_id}.{lane}",
                        samples=[readset.sample],
                        output_dependency=[output_json_path, output_html_path],
                        report_files=[output_json_path]
                    )
                )
                readset.report_files['fastp'] = [output_json_path]

                barcode_output_json_path = re.sub(".fastp.json", ".barcodes.fastp.json", output_json_path)
                barcode_output_html_path = re.sub(".fastp.html", ".barcodes.fastp.html", output_html_path)
                lane_jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                os.path.dirname(barcode_output_json_path),
                                remove=True
                            ),
                            fastp.basic_qc(
                                readset.index_fastq1,
                                readset.index_fastq2,
                                barcode_output_json_path,
                                barcode_output_html_path
                            )
                        ],
                        name=f"fastp.barcodes.{readset.name}.{self.run_id}.{lane}",
                        samples=[readset.sample],
                        output_dependency=[barcode_output_json_path, barcode_output_html_path],
                        report_files=[barcode_output_json_path]
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
                    name=f"blast." + readset.name + ".blast." + self.run_id + "." + lane,
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
                job.name = f"picard_mark_duplicates." + readset.name + ".dup." + self.run_id + "." + lane
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

    def check_sample_mixup(self):
        """
        Sample mixup analysis is performed to identify run processing data files from the same individual.
        This step can be used to identify doubles or compare samples in the same project id.
        Additionally, it can be used to compare current run processing data with an existing dataset.
        This enable us to find and compare samples from same individuals from two different runs.
        (ie. RNA and WGS data from two different projects)

        Currently it uses two commonly used tools (ie. NGSCheckmate and BAMixchecker) to identify sample mixups.

        NGSCheckmate only supports for human samples and FASTQ files (unaligned reads) are used (the tool supports BAM and VCF files.
        But the GenPipes only utilizes FASTQ files) as the input files.
        It supports various types of NGS data files including (but not limited to) whole genome sequencing (WGS),
        whole exome sequencing (WES), RNA-seq, single cell RNA-seq (scRNA) ChIP-seq, and targeted sequencing of various depths.
        NGScheckmate can identify sample mixup even from mixed data types (e.g. WES and RNA-seq, or RNA-seq and ChIP-seq).
        The analysis is performed for each project and all the lanes in a particular project
        Any type of mixed supported data can be used with in a project. NGSCheckMate uses depth-dependent
        correlation models of allele fractions of known single-nucleotide polymorphisms (SNPs) to identify samples from
        the same individual.
        For more information: https://github.com/parklab/NGSCheckMate
        The output files can be dound in sample_mixup_detection/NGSCheckMate/Homo_sapiens/PROJECT_ID

        BAMixchecker basically uses BAM files but BAM files need to be processed in-order to use as inputs to the tool.
        Depending on the sequencing type either 2 or no additional sub step will be executed automatically. Currently
        we have tested RNA-seq and DNA-seq data by ourselves. If the data are from RNA-seq
        below 2 steps will be done to pre-process data

        split_N_trim
        sambamba_merge_splitNtrim_files

        Pre-processed files will be then use as inputs to the BAMixchekcer. The output files can be found in
        sample_mixup_detection/BAMixChecker/projects/PROJECT_ID/SPECIES/BAMixChecker

        """

        jobs = []
        jobs.extend(self.checkmate_samplemixup_by_lane())
        jobs.extend(self.checkmate_samplemixup_by_run())
        jobs.extend(self.split_N_trim())
        jobs.extend(self.sambamba_merge_splitNtrim_files())
        jobs.extend(self.bamixchecker_samplemixup_by_lane())
        jobs.extend(self.bamixchecker_samplemixup_by_run())
        return jobs

    def checkmate_samplemixup_by_lane(self):
        """
        Check the sample mixup using NGSCheckmate by lane. FASTQ files are used as the inputs for the analysis
        Only human samples can be processed. If there are at least one human sample in the run, every sample in a lane
        is assume as human and run NGSCheckmate for all the samples in a lane
        If there is only one sample in a lane, it will be skipped
        """
        jobs = []

        #get all species in all the lanes and readsets and create a list

        #it doesn't check which species you have it always performs the analysis for humans. But there should be at
        # least one human sample in the whole run to perform the NGScheckmate analysis
        # as reference genome file is obtained from the event file

        # do not change below values. They are predefined by the tool
        output_file_names = ["output_all.txt", "output_corr_matrix.txt", "output_matched.txt"]

        all_ref_genome_paths = list(set([readset.reference_file for lane in self.lanes for readset in self.readsets[lane] if readset.reference_file]))
        genome = []
        if all_ref_genome_paths:
            genome = [fasta for fasta in all_ref_genome_paths if "Homo_sapiens" in fasta]

        # This job creates a ncm.config file used as an internal input to NGSCheckmate
        # do not change ncm_file name as the tool cannot identity any other name
        if len(genome) > 0:
            ncm_job = Job(
                [genome[0]],
                ["ncm.conf"],
                command="""echo -e "SAMTOOLS=samtools\nBCFTOOLS=bcftools\nREF={genome}" > {ncm_file}""".format(
                    genome = genome[0],
                    ncm_file="ncm.conf"
                )
            )
            ncm_job.name = "run_checkmate_ncf"
            jobs.append(ncm_job)

            # analysis is run for all the project ids separately. Since it only workd with human, doesn't loop over species.
            # and directly check whether the sample is from a human. No need to worry about the genome version since it
            # only takes the fastq but it checks with the defined reference genome in the ini or event file

            for lane in self.lanes:
                output_dir = "sample_mixup_detection/NGSCheckMate/Homo_sapiens/by_lane/"

                job_mkdir = bash.mkdir(os.path.join(output_dir, lane))
                samples = 0
                # create a list with output files. Output file names cannot be changed.
                output_files = [os.path.join(output_dir, lane) + "/" + sub for sub in output_file_names]

                filelist_list = []
                input_files = []
                lane_jobs = []

                filelist_path = os.path.join(output_dir, lane, "fastq_files_list_lane_" + lane +".txt")
                filelist_list.append(filelist_path)
                job_rm = bash.rm(filelist_path, force=True)
                job_touch = bash.touch(filelist_path)
                #filelist = None
                for readset in self.readsets[lane]:
                    samples += 1

                    #filelist = os.path.join(output_dir, lane, "fastq_files_list_lane_" + lane +".txt")

                    if (readset.fastq2):
                        filelist_job = Job(
                            [readset.fastq1, readset.fastq2],
                            [filelist_path],
                            command="""echo -e "{read1}\t{read2}\t{sample}" >> {file}""".format(
                                read1=readset.fastq1,
                                read2=readset.fastq2,
                                sample=readset.name,
                                file=filelist_path
                            )
                        )
                        input_files.append(readset.fastq1)
                        input_files.append(readset.fastq2)
                    else:
                        filelist_job = Job(
                            [readset.fastq1],
                            [filelist_path],
                            command="""echo -e "{read1}\t{sample}" >> {file}""".format(
                                read1=readset.fastq1,
                                sample=readset.name,
                                file=filelist_path
                            )
                        )
                        input_files.append(readset.fastq1)
                    filelist_job.samples = [readset.sample]
                    lane_jobs.append(filelist_job)
                job_filelist = concat_jobs(lane_jobs)

                if (samples > 1):
                    input_files.append("ncm.conf")
                    NGS_checkmate_job = ngscheckmate.fastq(
                        input_files,
                        output_files,
                        os.path.join(output_dir,lane),
                        filelist_path
                    )
                    job = concat_jobs(
                        [
                            job_mkdir, 
                            job_rm,
                            job_touch,
                            job_filelist,
                            NGS_checkmate_job
                        ]
                    )
                    job.name = "sample_mixup.ngscheckmate_by_lane_" + lane
                    jobs.append(job)
                else:
                    log.info("lane " + lane + " has only one sample... skipping NGSCheckmate analysis for this lane...")

        return self.throttle_jobs(jobs)

    def checkmate_samplemixup_by_run(self):
        """
        Check the sample mixup using NGSCheckmate. FASTQ files are used as the inputs for the analysis
        Only human samples can be processed.
        """
        jobs = []


        #get all species in all the lanes and readsets and create a list

        #it doesn't check which species you have it always performs the analysis for humans. But there should be at
        # least one human sample in the whole run to perform the NGScheckmate analysis
        # as reference genome file is obtained from the event file

        # do not change below values. They are predefined by the tool
        output_file_names = ["output_all.txt", "output_corr_matrix.txt", "output_matched.txt"]

        all_ref_genome_paths = list(set([readset.reference_file for lane in self.lanes for readset in self.readsets[lane] if readset.reference_file]))
        genome = [] 
        if all_ref_genome_paths:
            genome = [fasta for fasta in all_ref_genome_paths if "Homo_sapiens" in fasta]

        # This job creates a ncm.config file used as an internal input to NGSCheckmate
        # do not change ncm_file name as the tool cannot identity any other name
        if len(genome) > 0:

            # analysis is run for all the project ids separately. Since it only workd with human, doesn't loop over species.
            # and directly check whether the sample is from a human. No need to worry about the genome version since it
            # only takes the fastq but it checks with the defined reference genome in the ini or event file
            output_dir = "sample_mixup_detection/NGSCheckMate/Homo_sapiens/by_run/"

            job_mkdir = bash.mkdir(os.path.join(output_dir, self.run_id))
            samples = 0
            # create a list with output files. Output file names cannot be changed.

            filelist_list = []
            input_files = []
            run_jobs = []
            output_files = [os.path.join(output_dir, self.run_id) + "/" + sub for sub in output_file_names]
            filelist_path = os.path.join(output_dir, self.run_id, "fastq_files_list_run_" + self.run_id + ".txt")
            filelist_list.append(filelist_path)
            job_rm = bash.rm(filelist_path, force=True)
            job_touch = bash.touch(filelist_path)
            for lane in self.lanes:

                #filelist = None
                for readset in self.readsets[lane]:
                    samples += 1

                    #filelist = os.path.join(output_dir, lane, "fastq_files_list_lane_" + lane +".txt")

                    if (readset.fastq2):
                        filelist_job = Job(
                            [readset.fastq1, readset.fastq2],
                            [filelist_path],
                            command="""echo -e "{read1}\t{read2}\t{sample}" >> {file}""".format(
                                read1 = readset.fastq1,
                                read2 = readset.fastq2,
                                sample = readset.name,
                                file = filelist_path
                            )
                        )
                        input_files.append(readset.fastq1)
                        input_files.append(readset.fastq2)
                    else:
                        filelist_job = Job(
                            [readset.fastq1],
                            [filelist_path],
                            """echo -e "{read1}\t{sample}" >> {file}""".format(
                                read1=readset.fastq1,
                                sample=readset.name,
                                file=filelist_path
                            )
                        )
                        input_files.append(readset.fastq1)
                    filelist_job.samples = [readset.sample]
                    run_jobs.append(filelist_job)
                    job_filelist = concat_jobs(run_jobs)

            if(samples > 1):
                input_files.append("ncm.conf")
                NGS_checkmate_job = ngscheckmate.fastq(
                    input_files,
                    output_files,
                    os.path.join(output_dir,self.run_id),
                    filelist_path
                )

                job = concat_jobs(
                    [
                        job_mkdir,
                        job_rm,
                        job_touch,
                        job_filelist,
                        NGS_checkmate_job
                    ]
                )
                job.name = "sample_mixup.ngscheckmate_by_run_" + self.run_id
                jobs.append(job)
            else:
                log.info("Run " + self.run_id+ " has only one sample... skipping NGSCheckmate analysis for this run....")

        return self.throttle_jobs(jobs)

    def split_N_trim(self):
        """
        SplitNtrim. A [GATK](https://software.broadinstitute.org/gatk/) tool called SplitNCigarReads
        developed specially for RNAseq, which splits reads into exon segments (getting rid of Ns but
        maintaining grouping information) and hard-clip any sequences overhanging into the intronic regions.

        This analysis is only run for RNA-Seq data.
        """

        jobs = []

        nb_jobs = config.param('gatk_split_N_trim', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of haplotype jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for lane in self.lanes:
            lane_jobs = []
            for readset in [readset for readset in self.readsets[lane] if readset.bam]:
                sequence_dictionary = readset.dictionary_file

                if not (os.path.exists(sequence_dictionary)):
                    sequence_dictionary = config.param('DEFAULT', 'genome_dictionary', param_type='filepath', required=False)

                if (os.path.exists(sequence_dictionary) ):
                    bam = readset.bam + ".dup.bam"
                    if readset.is_rna:
                        # splitN trim only need to do for RNA
                        alignment_dir = os.path.join("sample_mixup_detection/BAMixChecker", "processed_BAMS",
                                                     readset.name)
                        split_dir = os.path.join(alignment_dir, "splitNtrim")
                        split_file_prefix = os.path.join(split_dir, readset.name + ".")

                        if nb_jobs == 1:
                            job = concat_jobs(
                                [
                                    bash.mkdir(alignment_dir),
                                    gatk.split_n_cigar_reads(
                                        bam,
                                        split_file_prefix + "sorted.dup.split.bam",
                                        fasta=readset.reference_file
                                    )
                                ],
                                name="gatk_split_N_trim." + readset.name,
                                samples=[readset.sample]
                            )
                            lane_jobs.append(job)

                        else:
                            unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(parse_sequence_dictionary_file(sequence_dictionary), nb_jobs - 1)
                            #  log.info(readset.name)
                            # Create one separate job for each of the first sequences
                            for idx, sequences in enumerate(unique_sequences_per_job):
                                # log.info(readset.name + "." + str(idx))
                                job = concat_jobs(
                                    [
                                        # Create output directory since it is not done by default by GATK tools
                                        bash.mkdir(split_dir),
                                        gatk.split_n_cigar_reads(
                                            bam,
                                            split_file_prefix + "sorted.dup.split." + str(idx) + ".bam",
                                            intervals=sequences,
                                            fasta=readset.reference_file
                                        )
                                    ],
                                    name="gatk_split_N_trim." + readset.name + "." + str(idx),
                                    samples=[readset.sample]
                                )
                                lane_jobs.append(job)

                            job = concat_jobs(
                                [
                                    bash.mkdir(alignment_dir),
                                    gatk.split_n_cigar_reads(
                                        bam,
                                        split_file_prefix + "sorted.dup.split.others.bam",
                                        exclude_intervals=unique_sequences_per_job_others,
                                        fasta=readset.reference_file
                                    )
                                ],
                                name="gatk_split_N_trim." + readset.name + ".others",
                                samples=[readset.sample],
                                removable_files=[split_file_prefix + "sorted.dup.split.others.bam"]
                            )
                            lane_jobs.append(job)
                else:
                    log.info("Please specify the Dictionary file path in an ini file. Skipping... ")
            jobs.extend(lane_jobs)
        return self.throttle_jobs(jobs)

    def sambamba_merge_splitNtrim_files(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [Sambamba] (http://lomereiter.github.io/sambamba/docs/sambamba-merge.html).
        """

        jobs = []

        nb_jobs = config.param('gatk_split_N_trim', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of haplotype jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for lane in self.lanes:
            lane_jobs =[]
            for readset in [readset for readset in self.readsets[lane] if readset.bam]:
                sequence_dictionary = readset.dictionary_file

                if not (os.path.exists(sequence_dictionary)):
                    sequence_dictionary = config.param('DEFAULT', 'genome_dictionary', param_type='filepath', required=False)

                if (os.path.exists(sequence_dictionary)):
                    if readset.is_rna:

                        alignment_directory = os.path.join("sample_mixup_detection/BAMixChecker", "processed_BAMS", readset.name)
                        split_file_prefix = os.path.join(alignment_directory, "splitNtrim", readset.name + ".")
                        output = os.path.join(alignment_directory, readset.name + ".sorted.dup.split.bam")

                        if nb_jobs > 1:
                            unique_sequences_per_job,_ = split_by_size(parse_sequence_dictionary_file(sequence_dictionary), nb_jobs - 1)

                            inputs = []
                            for idx,_ in enumerate(unique_sequences_per_job):
                                inputs.append(os.path.join(split_file_prefix + "sorted.dup.split." + str(idx) + ".bam"))
                            inputs.append(os.path.join(split_file_prefix + "sorted.dup.split.others.bam"))

                            job = sambamba.merge(
                                inputs,
                                output,
                                ini_section="sambamba_merge_splitNtrim_files"
                            )
                            job.name = "sambamba_merge_splitNtrim_files." + readset.name
                            job.samples = [readset.sample]
                            job.removable_files = [output]
                            lane_jobs.append(job)
                else:
                    log.info("Please specify the Dictionary file path in an ini file. Skipping... ")
            jobs.extend(lane_jobs)

        return self.throttle_jobs(jobs)

    def bamixchecker_samplemixup_by_lane(self):
        """
        Check the sample mixup using BAMixchecker. BAM files are used for the analysis
        """
        jobs = []

        species_list = list(set([readset.species for lane in self.lanes for readset in self.readsets[lane]]))

        for species in species_list:
            reference = None

            species_name = species.replace(" ", "_")
            species_name = species_name.replace(":", "_")
            species_name = re.sub('[^A-Za-z0-9_]+', '', species_name)
            for lane in self.lanes:
                samples = 0
                input_files = []
                output_dir = os.path.join("sample_mixup_detection/BAMixChecker",  species_name,  "by_lanes",lane)
                job_mkdir = bash.mkdir(output_dir)
                filelist = os.path.join(output_dir, "bam_files_list_" + species_name + "_" + lane + ".txt")
                job_rm = bash.rm(filelist, force=True)
                job_touch = bash.touch(filelist)
                lane_jobs = []

                for readset in [readset for readset in self.readsets[lane] if readset.bam]:
                    #and "Homo sapiens" in species

                    samples += 1
                    if readset.is_rna:
                        file_prefix = os.path.join("sample_mixup_detection/BAMixChecker", "processed_BAMS", readset.name, readset.name + ".sorted.dup.split.")
                        bam = file_prefix + "bam"
                    else:
                        file_prefix = os.path.join("sample_mixup_detection/BAMixChecker", "processed_BAMS", readset.name, readset.name + ".sorted.dup.")
                        bam = readset.bam + ".dup.bam"

                    filelist_job = Job(
                        [bam],
                        [filelist],
                        command="""echo -e "{bam}" >> {file}""".format(
                            bam=os.path.abspath(os.path.join(self.output_dir, bam)),
                            file=filelist
                        )
                    )
                    input_files.append(bam)
                    filelist_job.samples = [readset.sample]
                    lane_jobs.append(filelist_job)
                    job_filelist = concat_jobs(lane_jobs)

                    reference = readset.reference_file

                if samples > 1:
                    options = config.param('run_bamixchecker', 'options')
                    NonHumanSNPlist = config.param('run_bamixchecker', 'NonHumanSNPlist', required=False)

                    output_files = []
                    output_files.append(os.path.join(output_dir, "BAMixChecker", "Total_result.txt"))
                    output_files.append(os.path.join(output_dir, "BAMixChecker", "Matched_samples.txt"))

                    if "Homo_sapiens" in reference:
                        if "GRCh38" in reference or "hg38" in reference:
                            options = options + " -v hg38"
                        elif "GRCh37" in reference or "hg19" in reference:
                            options = options + " -v hg19"

                    elif "Mus_musculus" in reference:
                        if NonHumanSNPlist:
                            mouse_snp = NonHumanSNPlist
                            options = options + " --nhSNP " + mouse_snp

                    else:
                        if NonHumanSNPlist:
                            options = options + " --nhSNP " + NonHumanSNPlist

                    bammixchecker_job = bamixchecker.run_bam(
                        input_files,
                        output_files,
                        output_dir,
                        filelist,
                        reference,
                        options
                    )

                    if "Homo_sapiens" in reference or NonHumanSNPlist:
                        bammixchecker_job.samples = self.samples[lane]
                        job = concat_jobs(
                            [
                                job_mkdir,
                                job_rm,
                                job_touch,
                                job_filelist,
                                bammixchecker_job
                            ]
                        )
                        job.name = "sample_mixup.bamixchecker_by_lanes" + "_" + species_name + "_" + lane
                        jobs.append(job)

                else:
                    log.info(species_name +" of lane "+ lane + " has only one sample... skipping BAMixchekcer analysis for this lane....")

        return self.throttle_jobs(jobs)

    def bamixchecker_samplemixup_by_run(self):
        """
                Check the sample mixup using BAMixchecker. BAM files are used for the analysis
        """
        jobs = []

        species_list = list(set([readset.species for lane in self.lanes for readset in self.readsets[lane]]))

        for species in species_list:
            species_jobs = []
            reference = None

            species_name = species.replace(" ", "_")
            species_name = species_name.replace(":", "_")
            species_name = re.sub('[^A-Za-z0-9_]+', '', species_name)
            samples = 0
            input_files = []
            output_dir = os.path.join("sample_mixup_detection/BAMixChecker", species_name, "by_run", self.run_id)
            job_mkdir = bash.mkdir(output_dir)
            filelist = os.path.join(output_dir, "bam_files_list_" + species_name + "_" + self.run_id + ".txt")
            job_rm = bash.rm(filelist, force=True)
            job_touch = bash.touch(filelist)
            lane_jobs = []
            for lane in self.lanes:
                for readset in [readset for readset in self.readsets[lane] if readset.bam]:
                    #and "Homo sapiens" in species

                    samples += 1
                    if readset.is_rna:
                        file_prefix = os.path.join("sample_mixup_detection/BAMixChecker", "processed_BAMS", readset.name, readset.name + ".sorted.dup.split.")
                        bam = file_prefix + "bam"
                    else:
                        file_prefix = os.path.join("sample_mixup_detection/BAMixChecker", "processed_BAMS", readset.name, readset.name + ".sorted.dup.")
                        bam = readset.bam + ".dup.bam"

                    #bam = file_prefix + "recal.bam"

                    filelist_job = Job(
                        [bam],
                        [filelist],
                        command="""echo -e "{bam}" >> {file}""".format(
                            bam=os.path.abspath(os.path.join(self.output_dir, bam)),
                            file=filelist
                        )
                    )
                    input_files.append(bam)
                    filelist_job.samples = [readset.sample]
                    lane_jobs.append(filelist_job)
                    job_filelist = concat_jobs(lane_jobs)

                    reference = readset.reference_file

            if samples > 1:
                options = config.param('run_bamixchecker', 'options')
                NonHumanSNPlist = config.param('run_bamixchecker', 'NonHumanSNPlist', required=False)

                output_files = []
                output_files.append(os.path.join(output_dir, "BAMixChecker", "Total_result.txt"))
                output_files.append(os.path.join(output_dir, "BAMixChecker", "Matched_samples.txt"))

                if "Homo_sapiens" in reference:
                    if "GRCh38" in reference or "hg38" in reference:
                        options = options + " -v hg38"
                        bammixchecker_job = bamixchecker.run_bam(input_files, output_files, output_dir, filelist, reference, options)

                    elif "GRCh37" in reference or "hg19" in reference:
                        options = options + " -v hg19"
                        bammixchecker_job = bamixchecker.run_bam(input_files, output_files, output_dir, filelist, reference,
                                                                 options)
                # TO DO modify below code once prepare the mouse snp list
                elif "Mus_musculus" in reference:
                    if NonHumanSNPlist:
                        mouse_snp = NonHumanSNPlist
                        options = options + " --nhSNP " + mouse_snp
                        bammixchecker_job = bamixchecker.run_bam(input_files, output_files, output_dir, filelist,
                                                                 reference, options)
                else:
                    if NonHumanSNPlist:
                        options = options + " --nhSNP " + NonHumanSNPlist
                        bammixchecker_job = bamixchecker.run_bam(input_files, output_files, output_dir, filelist, reference,
                                                             options)

                if "Homo_sapiens" in reference or NonHumanSNPlist:
                    bammixchecker_job.samples = self.samples[lane]
                    job = concat_jobs([job_mkdir, job_rm, job_touch, job_filelist, bammixchecker_job])
                    job.name = "sample_mixup.bamixchecker_by_run" + "_" + species_name + "_" + self.run_id
                    jobs.append(job)

            else:
                log.info(species_name +" of run "+ self.run_id + " has only one sample... skipping BAMixchekcer analysis for this run....")

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

        return jobs
        # if self.args.type == 'mgit7':
        #     return jobs
        # else:
        #     return self.throttle_jobs(jobs)

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
                    name=f"md5." + readset.name + ".md5." + self.run_id + "." + lane,
                    samples=[readset.sample]
                )

                lane_jobs.append(job)

            if config.param('md5', 'one_job', required=False, param_type="boolean"):
                lane_job = concat_jobs(
                    lane_jobs,
                    name=f"md5.{self.run_id}.{lane}",
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

            self.generate_lane_json_report_file(lane)

            # metrics to JSON
            # Loop over all the steps of the pipeline
            step_list = [step for step in self.step_list if step.jobs]
            for step in step_list:
                report_step_jobs = []
                if step.name in ['basecall', 'fastq', 'index']:
                    step_report_files = []
                    for readset in self.readsets[lane]:
                        if step.name in readset.report_files:
                            step_report_files.extend(readset.report_files[step.name])
                    # With mgit7 platform running "slow-mode", set platform to mgig400
                    platform = "mgig400" if (self.args.type == "mgit7" and not self.args.splitbarcode_demux) else self.args.type
                    if step_report_files:
                        report_job = tools.run_processing_metrics_to_json(
                            self.run_validation_report_json[lane],
                            step.name,
                            platform,
                            list(set(step_report_files))
                        )
                        report_job.name = f"report." + step.name + "." + self.run_id + "." + lane
                        report_job.samples = self.samples[lane]
                        report_step_jobs.append(report_job)
                else:
                    for readset in self.readsets[lane]:
                        if step.name in readset.report_files:
                            # if step.name == 'fastp':
                            #     log.error(readset.report_files[step.name])
                            report_job = tools.run_processing_metrics_to_json(
                                self.run_validation_report_json[lane],
                                step.name,
                                self.args.type,
                                readset.report_files[step.name],
                                readset.name
                            )
                            report_job.name = f"report." + step.name + "." + readset.name + "." + self.run_id + "." + lane
                            report_job.samples = self.samples[lane]
                            # log.error(report_job.output_files)
                            report_step_jobs.append(report_job)

                lane_jobs.extend(self.throttle_jobs(report_step_jobs, f"{step.name}.{lane}"))

                # checkpoint file
                step_checkpoint_file = os.path.join(self.job_output_dir, f"checkpoint", step.name + "." + self.run_id + "." + lane + ".stepDone")
                step_checkpoint_job_dependencies = [output for job in step.jobs for output in job.output_files]
                if step_checkpoint_job_dependencies:
                    lane_jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(os.path.dirname(step_checkpoint_file)),
                                bash.touch(step_checkpoint_file)
                            ],
                            name=f"checkpoint." + step.name + "." + self.run_id + "." + lane,
                            input_dependency=step_checkpoint_job_dependencies,
                            output_dependency=[step_checkpoint_file],
                            samples=self.samples[lane]
                        )
                    )

            # loop over readsets and add file sizes of fastqs and bams to json
            for readset in self.readsets[lane]:

                if self.is_paired_end[lane]:
                    input_files = [readset.fastq1, readset.fastq2]
                else:
                    input_files = [readset.fastq1]
                if readset.bam:
                    input_files.extend([readset.bam + ".bam", readset.bam + ".bai"])

                size_job = tools.run_processing_file_sizes_to_json(
                    self.run_validation_report_json[lane],
                    input_files,
                    readset.name
                )
                size_job.name = f"report.file_sizes." + readset.name + "." + self.run_id + "." + lane
                size_job.samples = self.samples[lane]
                lane_jobs.append(size_job)

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
                os.path.basename(self.run_dir.rstrip('/')) + "-" + self.seqtype
            )

        jobs_to_concat.append(
            bash.mkdir(full_destination_folder)
        )

        jobs_to_concat.append(
            bash.cp(self.readset_file, full_destination_folder)
        )

        jobs_to_concat.append(
            bash.tar(
                [os.path.join(self.output_dir,
                              "job_output")],
                os.path.join(self.output_dir,
                             "job_output",
                             "job_output.tar.gz"),
                "".join(["-cpz --acls ",
                         "--exclude=\"",
                         os.path.join("job_output", "job_output.*\" "),
                         "--exclude=\"",
                         os.path.join("job_output", "copy", "copy.*.o\"")]),
                file_list=True,
                input_dependency=False)
        )

        jobs_to_concat.append(
            bash.cp(os.path.join(self.output_dir,
                                 "job_output",
                                 "job_output.tar.gz"),
                    os.path.join(full_destination_folder,
                                 "job_output")
                    )
        )

        jobs_to_concat.append(
            bash.cp(os.path.join(self.output_dir,
                                 "job_output",
                                 "job_output.list"),
                    os.path.join(full_destination_folder,
                                 "job_output")
                    )
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

            jobs_to_concat.append(
                concat_jobs(
                    [
                        Job(
                            command=copy_command_output_folder
                        ),
                        bash.touch(copy_output)
                    ],
                    input_dependency=inputs,
                    output_dependency=[copy_output]
                )
            )

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
            os.path.basename(self.run_dir.rstrip('/')) + "-" + self.seqtype
        )

        for lane in self.lanes:
            input = os.path.join(full_destination_folder, "copyCompleted." + lane + ".out")
            output = os.path.join(full_destination_folder, "notificationAssociation." + lane + ".out")

            notification_command = config.param('end_copy_notification', 'notification_command', required=False)
            if notification_command:
                job = Job(
                    [input],
                    [output],
                    name=f"end_copy_notification.{self.run_id}.{lane}"
                )
                job.command = notification_command.format(
                    technology=config.param('end_copy_notification', 'technology'),
                    output_dir=self.output_dir,
                    run_name=os.path.basename(self.run_dir.rstrip('/')),
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
                    os.path.basename(self.run_dir.rstrip('/')) + "-" + self.seqtype
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
                "input_files" : [os.path.relpath(input_file, self.output_dirs[lane]["report_directory"]) for input_file in job.input_files],
                "output_files" : [os.path.relpath(output_file, self.output_dirs[lane]["report_directory"]) for output_file in job.output_files],
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
            return str(int(self.json_flag_hash[lane]['Read1']) - 1)

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
            return str(int(self.json_flag_hash[lane]['Read2']) - 1)
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
        Returns the number of cycles for each index reads of the run.
        """
        for read in [read for read in self.read_infos if read.is_index]:
            if read.number in [1, 2]:
                index1cycles = read.nb_cycles
            if read.number in [3, 4]:
                index2cycles = read.nb_cycles
        return [index1cycles, index2cycles]

    def get_indexreverse(self, lane):
        """
        Returns the reverse complement status for each index reads of the run
        if available in the current platform.
        """
        if self.args.type == "illumina":
            [index1reverse, index2reverse] = self.get_illumina_indexreverse(lane)
            return [index1reverse, index2reverse]
        else:
            _raise(SanitycheckError(
                "Unsupported reverse index from RunInfo.xml on platform: " + self.args.type))

    def get_illumina_indexreverse(self, lane):
        """
        Returns the reverse complement status for each index reads of the run.
        """
        for read in [read for read in self.read_infos if read.is_index]:
            if read.number in [1, 2]:
                index1reverse = read.is_reverse_complement
            if read.number in [3, 4]:
                index2reverse = read.is_reverse_complement
        return [index1reverse, index2reverse]

    def get_sequencer_index_length(self):
        """
        Returns the total number of index cycles of the run.
        """
        return sum(index_read.nb_cycles for index_read in [read for read in self.read_infos if read.is_index])

    def get_sequencer_minimum_read_length(self, lane):
        """
        Returns the minimum number of cycles of a real read (not indexed).
        """
        return min([self.read1cycles[lane], self.read2cycles[lane]])

    # TODO Obsolete... Is not called in this file. Must make sure it is not
    # imported and used somewhere else before deleting the function.
    def has_single_index(self, lane):
        """
        Returns True when there is at least one sample on the lane that doesn't use double-indexing or we only have
        one read of indexes.

        Obsolete
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

        all_indexes = []
        for readset in self.readsets[lane]:
            all_indexes += readset.indexes

        # loop on all index reads, to compare with samples index length
        for i in range(len(run_index_lengths)):
            min_sample_index_length = 0
            try:
                min_sample_index_length = min(len(index['INDEX'+str(i+1)])
                                              for index in all_indexes
                                              if (len(index) > i and len(index['INDEX'+str(i+1)])))
            except ValueError:
                pass  # we don't have a sample with this Ith index read, use the 0 already set

            empty_index_list = [ index for index in all_indexes if (len(index) <= i or len(index['INDEX'+str(i+1)]) == 0) ]

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
        """
        Validate all index sequences against each other to ensure they aren't
        in collision according to the chosen number of mismatches parameter.
        """
        min_allowed_distance = (2 * self.number_of_mismatches) + 1
        index_lengths = self.get_smallest_index_length(lane)

        validated_indexes = []
        collisions = []

        for readset in self.readsets[lane]:
            for current_index in readset.indexes:
                if int(self.index2cycles[lane]):
                    if 'mgi' in self.args.type:
                        current_barcode = current_index['INDEX2']+current_index['INDEX1']
                    else:
                        current_barcode = current_index['INDEX1']+current_index['INDEX2']
                else:
                    if current_index['INDEX1']:
                        current_barcode = current_index['INDEX1']
                    else:
                        current_barcode = current_index['INDEX2']
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
        if int(self.index2cycles[lane]):
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
        This Flag file is used to call client_linux (i.e wrapper of splitBarcode and more)
        """
        json_flag_file = self.json_flag_files[lane]

        # get the barcode names & sequences to add in the JSON flag file
        all_indexes = {}
        for readset in self.readsets[lane]:
            for index_dict in readset.index:
                all_indexes[readset.index_name] = index_dict
        with open(json_flag_file, 'r') as json_fh:
            json_flag_content = json.load(json_fh)

        if self.is_dual_index[lane]:
            json_flag_content['speciesBarcodes'] = dict([(index_name, index_dict['INDEX2']+index_dict['INDEX1']) for index_name, index_dict in all_indexes.items()])
        else:
            json_flag_content['speciesBarcodes'] = dict([(index_name, index_dict['INDEX1']) for index_name, index_dict in all_indexes.items()])
        json_flag_content['barcodeStartPos'] = json_flag_content['TotalCycle'] - json_flag_content['barcodeLength'] + 1
        json_flag_content['SpeciesMismatch'] = self.number_of_mismatches

        with open(json_flag_file, 'w') as out_json_fh:
            json.dump(json_flag_content, out_json_fh, indent=4)

        log.info("BARCODES added in FLAG file : " + json_flag_file)

    def create_barcode_file(self, lane):
        """
        create the so-called barcode file used to call MGI splitBarcode directly
        --> could be useless in the end because it seems splitBarcode could also use the flag file (still need to be tested)
        """

        barcode_file = self.barcode_files[lane]

        # get the barcode names & sequences to add in the JSON flag file
        all_indexes = {}
        for readset in self.readsets[lane]:
            for index_dict in readset.index:
                all_indexes[readset.index_name] = index_dict

        if self.is_dual_index[lane]:
            barcodes = dict([(index_name, index_dict['INDEX2']+"\t"+index_dict['INDEX1']) for index_name, index_dict in all_indexes.items()])
        else:
            barcodes = dict([(index_name, index_dict['INDEX1']) for index_name, index_dict in all_indexes.items()])

        with open(barcode_file, 'w') as barcode_fh:
            for barcode in barcodes.keys():
                barcode_fh.write(barcode + "\t" + barcodes[barcode] + "\n")

        log.info("BARCODE FILE created : " + barcode_file)

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

         # barcode validation
        if re.search("I", mask) and not self.args.allow_barcode_collision:
            self.validate_barcodes(lane)

        overmask = ""
        overindex1 = None
        overindex2 = None

        self._umi = False

        # IDT - UMI9 in index2
        if re.search(",I17", mask):
            self._umi = True
            overmask = re.sub(",I17,", ",I8n*,", mask)
            overindex1=8
            overindex2=8

        maxindex1cycles = min(int(self.index1cycles[lane]), max([len(index['INDEX1']) for readset in self.readsets[lane] for index in readset.index]))
        maxindex2cycles = min(int(self.index2cycles[lane]), max([len(index['INDEX2']) for readset in self.readsets[lane] for index in readset.index]))
        pad1 = "n*" if maxindex1cycles != self.index1cycles[lane] else ""
        pad2 = "n*" if maxindex2cycles != self.index2cycles[lane] else ""

        split_mask = mask.split(",")

        # If SINGLEINDEX only
        if "DUALINDEX" not in set([readset.index_type for readset in self.readsets[lane]]):
            if re.search("I", mask.split(",")[1]) and len([readset for readset in self.readsets[lane] if readset.library_type in ["tenX_sc_RNA_v1", "TELL-Seq", "SHARE-Seq_ATAC", "SHARE-Seq_RNA"]]):
                # R2 is I1
                overmask=','.join([split_mask[0], 'Y'+self.index1cycles[lane], 'I'+str(maxindex2cycles)+pad2, split_mask[3]])
                overindex1=0
                overindex2=maxindex2cycles
                self._bcl2fastq_extra_option="--mask-short-adapter-reads 8"
            elif re.search("I", mask.split(",")[2]):
                # R2 is I2
                overmask=','.join([split_mask[0], 'I'+str(maxindex1cycles)+pad1, 'Y'+self.index2cycles[lane], split_mask[3]])
                overindex1=maxindex1cycles
                overindex2=0
                self._bcl2fastq_extra_option="--mask-short-adapter-reads 8"
            else:
                # read 3 is not index
                overmask=','.join([split_mask[0], 'I'+str(maxindex1cycles)+pad1, split_mask[2]])
                overindex1=maxindex1cycles
                overindex2=0

        # If DUALINDEX or MIX
        else:
            if re.search("I", mask.split(",")[2]):
                # read 3 is sequenced as index
                overmask=','.join([split_mask[0], 'I'+str(maxindex1cycles)+pad1, 'I'+str(maxindex2cycles)+pad2, split_mask[3]])
                overindex1=maxindex1cycles
                overindex2=maxindex2cycles
            else:
                # read 3 is not index
                overmask=','.join([split_mask[0], 'I'+str(maxindex1cycles)+pad1, split_mask[2]])
                overindex1=maxindex1cycles
                overindex2=0

        final_mask = mask if overmask == "" else overmask
        final_index1 = self.index1cycles[lane] if overindex1 is None else overindex1
        final_index2 = self.index2cycles[lane] if overindex2 is None else overindex2

        self._mask[lane] = config.param('fastq_illumina', 'overmask') if config.param('fastq_illumina', 'overmask', required=False, param_type='string') else final_mask
        self._index1cycles[lane] = config.param('fastq_illumina', 'overindex1') if config.param('fastq_illumina', 'overindex1', required=False, param_type='int') else final_index1
        self._index2cycles[lane] = config.param('fastq_illumina', 'overindex2') if config.param('fastq_illumina', 'overindex2', required=False, param_type='int') else final_index2

        index_lengths = self.get_smallest_index_length(lane)
        for readset in sorted(self.readsets[lane], key=lambda item: int(item.sample_number)):

            for idx, readset_index in enumerate(readset.indexes):
                # Barcode sequence should only match with the barcode cycles defined in the mask
                # so we adjust the length of the index sequences accordingly for the "Sample_Barcode" field
                if int(self.index2cycles[lane]):
                    sample_barcode = readset_index['INDEX1'][0:index_lengths[0]] + readset_index['INDEX2'][0:index_lengths[1]]
                else:
                    sample_barcode = readset_index['INDEX1'][0:index_lengths[0]]
                if self.last_index < len(sample_barcode):
                    sample_barcode = sample_barcode[0:self.last_index]
                if self.first_index > 1:
                    sample_barcode = sample_barcode[self.first_index-1:]

                # Reverse complement the indexes that are marked in RunInfo v6
                # Only relevent to the NovaseqX so far
                if self.index1reverse[lane]:
                    readset_index['INDEX1'] = readset_index['INDEX1'].replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c").upper()[::-1]
                if self.index2reverse[lane]:
                    readset_index['INDEX2'] = readset_index['INDEX2'].replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c").upper()[::-1]

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
                    "Description": readset.index_name + ' - ' + readset.protocol + ' - ' + readset.library_source,
                    "Control": readset.control,
                    "Recipe": readset.recipe,
                    "Operator": readset.operator,
                    "Sample_Project": "Project_" + readset.project_id
                }
                writer.writerow(csv_dict)

        if self.umi or len(self.readsets[lane]) == 1:
            casava_sample_sheet_noindex = re.sub(".indexed.", ".noindex.", csv_file)
            writer = csv.DictWriter(
                open(casava_sample_sheet_noindex, 'w'),
                delimiter=str(','),
                fieldnames=csv_headers
            )
            writer.writerow(
                {
                    "FCID": self.readsets[lane][0].flow_cell,
                    "Lane": lane,
                    "Sample_ID": "Sample_ALL" if self.umi else "Sample_"+self.readsets[lane][0].name,
                    "Sample_Name": "ALL" if self.umi else self.readsets[lane][0].name,
                    "SampleRef": "",
                    "Index": "",
                    "Index2": "",
                    "Description": "" if self.umi else self.readsets[lane][0].index_name + ' - ' + self.readsets[lane][0].protocol + ' - ' + self.readsets[lane][0].library_source,
                    "Control": "N",
                    "Recipe": "",
                    "Operator": "",
                    "Sample_Project": "Project_ALL" if self.umi else "Project_" + self.readsets[lane][0].project_id
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
                # so we adjust the length of the index sequences accordingly for the "Sample_Barcode" field
                if int(self.index2cycles[lane]):
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
                        "barcode": readset.index_name,
                        "barcode_sequence": ','.join([readset_index['BARCODE_SEQUENCE'] for readset_index in readset.indexes]),
                        "pct_on_index_in_lane": None,
                        "pct_of_the_lane": None,
                        "pct_perfect_barcode": None,
                        "pct_one_mismatch_barcode": None,
                        "pf_clusters": None
                    },
                    "qc": {
                        "yield": None,
                        "pct_q30_bases": None,
                        "avg_qual": None,
                        "duplicate_rate": None,
                        "nb_reads": None
                    },
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
                        "adapter_dimers": None,
                        "mean_coverage": None,
                        "aligned_dup_rate": None,
                        "inferred_sex": None,
                        "sex_concordance": None
                    }
                }
            )

        # Adding the MultiQC input file list
        # self.report_hash[lane]["multiqc_inputs"] = []
        step_list = [step for step in self.step_list if step.jobs]
        self.report_hash[lane]["multiqc_inputs"] = list(set([report_file for step in step_list for job in step.jobs for report_file in job.report_files if f"ligned.{lane}" in report_file]))
        self.report_hash[lane]["multiqc_inputs"].append(os.path.join(self.output_dirs[lane]["report_directory"], f"{self.run_id}.{lane}.run_validation_report.json"))
        self.report_hash[lane]["metrics_report_url"] = f"https://datahub-297-p25.p.genap.ca/Freezeman_validation/{self.year}/{self.run_id}.report.html"

        if not os.path.exists(os.path.join(self.output_dir, os.path.dirname(self.run_validation_report_json[lane]))):
            os.makedirs(os.path.join(self.output_dir, os.path.dirname(self.run_validation_report_json[lane])))
        if not os.path.exists(os.path.join(self.output_dir, self.run_validation_report_json[lane])) or self.force_jobs:
            with open(os.path.join(self.output_dir, self.run_validation_report_json[lane]), 'w') as out_json:
                json.dump(self.report_hash[lane], out_json, indent=4)

    def generate_basecall_outputs(self, lane):
        basecall_outputs = []
        postprocessing_jobs = []
        unaligned_dir = self.output_dirs[lane][f"Unaligned.{lane}_directory"]
        basecall_dir = os.path.join(unaligned_dir, "basecall")

        index_lengths = self.get_smallest_index_length(lane)

        for readset in self.readsets[lane]:
            readset_r1_outputs = []
            readset_r2_outputs = []
            for index in readset.indexes:
                readset_r1_outputs.extend(
                    [
                        os.path.join(basecall_dir, self.run_id, f"L0{lane}", self.flowcell_id +  "_L0" + lane + "_" + index['INDEX_NAME'] + "_1.fq.gz"),
                        os.path.join(basecall_dir, self.run_id, f"L0{lane}", self.flowcell_id +  "_L0" + lane + "_" + index['INDEX_NAME'] + "_1.fq.fqStat.txt")
                    ]
                )
                if readset.run_type == "PAIRED_END":
                    readset_r2_outputs.extend(
                        [
                            os.path.join(basecall_dir, self.run_id, f"L0{lane}", self.flowcell_id +  "_L0" + lane + "_" + index['INDEX_NAME'] + "_2.fq.gz"),
                            os.path.join(basecall_dir, self.run_id, f"L0{lane}", self.flowcell_id +  "_L0" + lane + "_" + index['INDEX_NAME'] + "_2.fq.fqStat.txt")
                        ]
                    )

                # replace with metrics job that produces per-lane output, instead of per-sample output
                #summary_report = re.sub("_1.fq.gz", "_summaryReport.html", readset_r1_outputs[0])

                #summary_report_job = run_processing_tools.mgi_summary_report(
                #    readset_r1_outputs[1],
                #    readset_r2_outputs[1] if readset.run_type == "PAIRED_END" else None,
                #    summary_report,
                #    index['INDEX_NAME']
                #    )
                #summary_report_job.name = "mgi_summary_report." + index['INDEX_NAME']
                #summary_report_job.samples=self.samples[lane]
                #postprocessing_jobs.append(summary_report_job)

            # If True, then merge the 'Undetermined' reads
            if self.merge_undetermined[lane]:
                readset_r1_outputs.extend(
                    [
                        os.path.join(basecall_dir, self.run_id, f"L0{lane}", self.flowcell_id +  "_L0" + lane + "_undecoded_1.fq.gz"),
                        os.path.join(basecall_dir, self.run_id, f"L0{lane}", self.flowcell_id +  "_L0" + lane + "_undecoded_1.fq.fqStat.txt")
                    ]
                )
                if readset.run_type == "PAIRED_END":
                    readset_r2_outputs.extend(
                        [
                            os.path.join(basecall_dir, self.run_id, f"L0{lane}", self.flowcell_id +  "_L0" + lane + "_undecoded_2.fq.gz"),
                            os.path.join(basecall_dir, self.run_id, f"L0{lane}", self.flowcell_id +  "_L0" + lane + "_undecoded_2.fq.fqStat.txt")
                        ]
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
                                    self.awk_read_1_processing_command(lane, readset)
                                ),
                                bash.gzip(
                                    None,
                                    readset.fastq1
                                )
                            ]
                        ),
                        samtools.quickcheck(
                            readset.fastq1,
                            None,
                            "-u"
                            )
                    ],
                    name=f"fastq_convert.R1." + readset.name + "." + self.run_id + "." + lane,
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
                            ),
                            samtools.quickcheck(
                                readset.fastq2,
                                None,
                                "-u"
                                )
                        ],
                        input_dependency=readset_r2_outputs,
                        output_dependency=outputs,
                        name=f"fastq_convert.R2." + readset.name + "." + self.run_id + "." + lane,
                        samples=self.samples[lane]
                    )
                )
        # Produce summaryReport for lane
        lane_basecall_dir = os.path.join(basecall_dir, self.run_id, f"L0{lane}")
        metrics_dir = os.path.join(self.run_dir, f"L0{lane}", "metrics") 
        
        if self.is_paired_end[lane]:
            PE = "-p"
        else:
            PE = None

        postprocessing_jobs.append(
            concat_jobs(
                [
                    run_processing_tools.mgi_lane_summary_report(
                        metrics_dir,
                        lane_basecall_dir,
                        self.run_id,
                        PE
                        )
                ],
                name="mgi_lane_summary_report." + f"L0{lane}",
                input_dependency=basecall_outputs
                )
            )
        # Process undetermined reads fastq files
        unmatched_R1_fastq = os.path.join(basecall_dir, self.run_id, f"L0{lane}", self.flowcell_id +  "_L0" + lane + "_undecoded_1.fq.gz")
        if unmatched_R1_fastq not in basecall_outputs:
            basecall_outputs.append(unmatched_R1_fastq)
        postprocessing_jobs.append(
            concat_jobs(
                [
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
                    ),
                    samtools.quickcheck(
                        os.path.join(unaligned_dir, "Undetermined_S0_L00" + lane + "_R1_001.fastq.gz"),
                        None,
                        "-u"
                        )
                ],
                name=f"fastq_convert.R1.unmatched.{self.run_id}.{lane}",
                samples=self.samples[lane]
            )
        )
        
        if self.is_paired_end[lane]:
            unmatched_R2_fastq = os.path.join(basecall_dir, self.run_id, f"L0{lane}", self.flowcell_id +  "_L0" + lane + "_undecoded_2.fq.gz")
            unaligned_i1 = os.path.join(unaligned_dir, "Undetermined_S0_L00" + lane + "_I1_001.fastq.gz")
            unexpected_barcode_counts_i1 = re.sub(".fastq.gz", ".counts.txt", unaligned_i1)
            if unmatched_R2_fastq not in basecall_outputs:
                basecall_outputs.append(unmatched_R2_fastq)
            outputs = [
                os.path.join(unaligned_dir, "Undetermined_S0_L00" + lane + "_R2_001.fastq.gz"),
                unaligned_i1
            ]
            if self.is_dual_index[lane]:
                unaligned_i2 = os.path.join(unaligned_dir, "Undetermined_S0_L00" + lane + "_I2_001.fastq.gz")
                unexpected_barcode_counts_i2 = re.sub(".fastq.gz", ".counts.txt", unaligned_i2)
                outputs.append(unaligned_i2)
            else:
                unaligned_i2 = False
            postprocessing_jobs.append(
                concat_jobs(
                    [
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
                        ),
                        samtools.quickcheck(
                            os.path.join(unaligned_dir, "Undetermined_S0_L00" + lane + "_R2_001.fastq.gz"),
                            None,
                            "-u"
                        )
                    ],
                    input_dependency=[unmatched_R2_fastq],
                    output_dependency=outputs,
                    name=f"fastq_convert.R2.unmatched.{self.run_id}.{lane}",
                    samples=self.samples[lane]
                )
            )

            if unaligned_i1:
                index1length = index_lengths[1] if unaligned_i2 else index_lengths[0]
                if index1length:
                    postprocessing_jobs.append(
                        pipe_jobs(
                            [
                                bash.cat(
                                    unaligned_i1,
                                    None,
                                    zip=True
                                ),
                                bash.awk(
                                    None,
                                    unexpected_barcode_counts_i1,
                                    f"'NR%4==2 {{print substr($0,1,{index1length})}}' | sort | uniq -c | sort -nr"
                                )
                            ],
                            name="fastq_countbarcodes.I1.unmatched." + self.run_id + "." + lane,
                            samples=self.samples[lane]
                        )
                    )
            if unaligned_i2:
                index1length = index_lengths[1]
                index2length = index_lengths[0]
                if index2length:
                    postprocessing_jobs.append(
                        pipe_jobs(
                            [
                                bash.cat(
                                    unaligned_i2,
                                    None,
                                    zip=True
                                ),
                                bash.awk(
                                    None,
                                    unexpected_barcode_counts_i2,
                                    f"'NR%4==2 {{print substr($0,1,{index2length})}}' | sort | uniq -c | sort -nr"
                                )
                            ],
                            name="fastq_countbarcodes.I2.unmatched." + self.run_id + "." + lane,
                            samples=self.samples[lane]
                        )
                    )
                unexpected_barcode_counts_i1i2 = re.sub("_I1_", "_I1I2_", unexpected_barcode_counts_i1)
                postprocessing_jobs.append(
                    concat_jobs(
                        [
                            bash.paste(
                                None,
                                unexpected_barcode_counts_i1i2,
                                options=f"-d '' <(zcat {unaligned_i1} | awk 'NR%4==2 {{print substr($0,1,{index1length})}}') <(zcat {unaligned_i2} | awk 'NR%4==2 {{print substr($0,1,{index2length})}}') | sort | uniq -c | sort -nr"
                            )
                        ],
                        input_dependency=[
                            unaligned_i1,
                            unaligned_i2
                            ],
                        name="fastq_countbarcodes.I1I2.unmatched." + self.run_id + "." + lane,
                        samples=self.samples[lane]
                    )
                )
    
            if unaligned_i2:
                unexpected_barcode_counts = unexpected_barcode_counts_i1i2
            else:
                unexpected_barcode_counts = unexpected_barcode_counts_i1
            unexpected_barcode_matches = re.sub(".counts.txt", ".match_table.tsv", unexpected_barcode_counts)
    
            job = run_processing_tools.match_undetermined_barcodes(
                    unexpected_barcode_counts,
                    unexpected_barcode_matches
                    )
            job.name = 'fastq_match_undetermined_barcodes.' + self.run_id + "." + lane
            job.samples = self.samples[lane]
            postprocessing_jobs.append(job)
        
        return basecall_outputs, postprocessing_jobs

    def generate_bcl2fastq_outputs(self, lane):
        bcl2fastq_outputs = []
        postprocessing_jobs = []
        cleanjob_deps = []
        
        output_dir = self.output_dirs[lane][f"Unaligned.{lane}_directory_tmp"]

        bcl2fastq_outputs.extend(
            [
                os.path.join(output_dir, "Reports"),
                os.path.join(output_dir, "Stats"),
                os.path.join(output_dir, f"Undetermined_S0_L00{lane}_R1_001.fastq.gz"),
                os.path.join(output_dir, f"Undetermined_S0_L00{lane}_I1_001.fastq.gz")
            ]
        )
        if self.is_paired_end[lane]:
            bcl2fastq_outputs.append(os.path.join(output_dir, f"Undetermined_S0_L00{lane}_R2_001.fastq.gz"))
        if self.is_dual_index[lane]:
            bcl2fastq_outputs.append(os.path.join(output_dir, f"Undetermined_S0_L00{lane}_I2_001.fastq.gz"))

        sample_count = 0
        for readset in self.readsets[lane]:
            readset_r1_outputs = []
            readset_r2_outputs = []
            readset_i1_outputs = []
            readset_i2_outputs = []
            readset_r3_outputs = [] # For Haloplex-like libraries

            for index in readset.indexes:
                sample_count += 1
                readset_index_outdir = os.path.join(output_dir, f"Project_{readset.project_id}", f"Sample_{index['SAMPLESHEET_NAME']}")
                readset_r1_outputs.append(
                    os.path.join(readset_index_outdir, f"{index['SAMPLESHEET_NAME']}_S{sample_count}_L00{lane}_R1_001.fastq.gz")
                )
                readset_i1_outputs.append(
                    os.path.join(readset_index_outdir, f"{index['SAMPLESHEET_NAME']}_S{sample_count}_L00{lane}_I1_001.fastq.gz")
                )
                if self.is_paired_end[lane]:
                    readset_r2_outputs.append(
                        os.path.join(readset_index_outdir, f"{index['SAMPLESHEET_NAME']}_S{sample_count}_L00{lane}_R2_001.fastq.gz")
                    )
                if self.is_dual_index[lane]:
                    readset_i2_outputs.append(
                        os.path.join(readset_index_outdir, f"{index['SAMPLESHEET_NAME']}_S{sample_count}_L00{lane}_I2_001.fastq.gz")
                    )
                # There is an R3 for all Haloplex-like libraries
                if ''.join(i for i in self.mask if not i.isdigit()) in ["Y,I,Y,Y", "Y,Y,I,Y"]:
                    readset_r3_outputs.append(
                        os.path.join(readset_index_outdir, f"{index['SAMPLESHEET_NAME']}_S{sample_count}_L00{lane}_R3_001.fastq.gz")
                    )
                

            # If True, then merge the 'Undetermined' reads
            if self.merge_undetermined[lane]:
                readset_r1_outputs.append(
                    os.path.join(output_dir, f"Undetermined_S0_L00{lane}_R1_001.fastq.gz")
                )
                readset_i1_outputs.append(
                    os.path.join(output_dir, f"Undetermined_S0_L00{lane}_I1_001.fastq.gz")
                )
                if self.is_paired_end[lane]:
                    readset_r2_outputs.append(
                        os.path.join(output_dir, f"Undetermined_S0_L00{lane}_R2_001.fastq.gz")
                    )
                if self.is_dual_index[lane]:
                    readset_i2_outputs.append(
                        os.path.join(output_dir, f"Undetermined_S0_L00{lane}_I2_001.fastq.gz")
                    )
                # R3 for Haloplex-like libraries
                if ''.join(i for i in self.mask if not i.isdigit()) in ["Y,I,Y,Y", "Y,Y,I,Y"]:
                    readset_r3_outputs.append(
                        os.path.join(output_dir, f"Undetermined_S0_L00{lane}_R3_001.fastq.gz")
                    )

            bcl2fastq_outputs.extend(readset_r1_outputs)
            if not len(readset_r1_outputs) == 1:
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
                                    bash.gzip(
                                        None,
                                        readset.fastq1
                                    )
                                ]
                            )
                        ],
                        name=f"fastq_concat.R1." + readset.name + "." + self.run_id + "." + lane,
                        samples=[readset.sample]
                    )
                )
            else:
                postprocessing_jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(os.path.dirname(readset.fastq1)),
                            bash.mv(readset_r1_outputs[0], readset.fastq1)
                        ],
                        name=f"fastq_move.R1." + readset.name + "." + self.run_id + "." + lane,
                        samples=[readset.sample]
                    )
                )
            cleanjob_deps.append(readset.fastq1)

            bcl2fastq_outputs.extend(readset_i1_outputs)
            if not len(readset_i1_outputs) == 1:
                postprocessing_jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(os.path.dirname(readset.index_fastq1)),
                            pipe_jobs(
                                [
                                    bash.cat(
                                        readset_i1_outputs,
                                        None,
                                        zip=True
                                    ),
                                    bash.gzip(
                                        None,
                                        readset.index_fastq1
                                    )
                                ]
                            )
                        ],
                        name=f"fastq_concat.I1." + readset.name + "." + self.run_id + "." + lane,
                        samples=[readset.sample]
                    )
                )
            else:
                postprocessing_jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(os.path.dirname(readset.index_fastq1)),
                            bash.mv(readset_i1_outputs[0], readset.index_fastq1)
                        ],
                        name=f"fastq_move.I1." + readset.name + "." + self.run_id + "." + lane,
                        samples=[readset.sample]
                    )
                )
            cleanjob_deps.append(readset.index_fastq1)

            if self.is_paired_end[lane]:
                bcl2fastq_outputs.extend(readset_r2_outputs)
                if not len(readset_r2_outputs) == 1:
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
                                        bash.gzip(
                                            None,
                                            readset.fastq2
                                        )
                                    ]
                                )
                            ],
                            name=f"fastq_concat.R2." + readset.name + "." + self.run_id + "." + lane,
                            samples=[readset.sample]
                        )
                    )
                else:
                    postprocessing_jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(os.path.dirname(readset.fastq2)),
                                bash.mv(readset_r2_outputs[0], readset.fastq2)
                            ],
                            name=f"fastq_move.R2." + readset.name + "." + self.run_id + "." + lane,
                            samples=[readset.sample]
                        )
                    )
                cleanjob_deps.append(readset.fastq2)

            if self.is_dual_index[lane]:
                bcl2fastq_outputs.extend(readset_i2_outputs)
                if not len(readset_i2_outputs) == 1:
                    postprocessing_jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(os.path.dirname(readset.index_fastq2)),
                                pipe_jobs(
                                    [
                                        bash.cat(
                                            readset_i2_outputs,
                                            None,
                                            zip=True
                                        ),
                                        bash.gzip(
                                            None,
                                            readset.index_fastq2
                                        )
                                    ]
                                )
                            ],
                            name="fastq_concat.I2." + readset.name + "." + self.run_id + "." + lane,
                            samples=[readset.sample]
                        )
                    )
                else:
                    postprocessing_jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(os.path.dirname(readset.index_fastq2)),
                                bash.mv(readset_i2_outputs[0], readset.index_fastq2)
                            ],
                            name="fastq_move.I2." + readset.name + "." + self.run_id + "." + lane,
                            samples=[readset.sample]
                        )
                    )
                cleanjob_deps.append(readset.index_fastq2)

            # R3 for Haloplex-like libraries
            if ''.join(i for i in self.mask if not i.isdigit()) in ["Y,I,Y,Y", "Y,Y,I,Y"]:
                bcl2fastq_outputs.extend(readset_r3_outputs)
                if not len(readset_r3_outputs) == 1:
                    postprocessing_jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(os.path.dirname(re.sub("_R1_", "_R3_", readset.fastq1))),
                                pipe_jobs(
                                    [
                                        bash.cat(
                                            readset_r3_outputs,
                                            None,
                                            zip=True
                                        ),
                                        bash.gzip(
                                            None,
                                            re.sub("_R1_", "_R3_", readset.fastq1)
                                        )
                                    ]
                                )
                            ],
                            name="fastq_concat.R3." + readset.name + "." + self.run_id + "." + lane,
                            samples=[readset.sample]
                        )
                    )
                else:
                    postprocessing_jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(os.path.dirname(re.sub("_R1_", "_R3_", readset.fastq1))),
                                bash.mv(readset_r3_outputs[0], re.sub("_R1_", "_R3_", readset.fastq1))
                            ],
                            name="fastq_move.R3." + readset.name + "." + self.run_id + "." + lane,
                            samples=[readset.sample]
                        )
                    )

            # For HaloPlex-like libraries, do the necessary name swapping
            if ''.join(i for i in self.mask if not i.isdigit()) == "Y,I,Y,Y":
                postprocessing_jobs.append(
                    concat_jobs(
                        [
                            bash.mv(
                                readset.fastq2,
                                readset.index_fastq2
                            ),
                            bash.mv(
                                re.sub("_R1_", "_R3_", readset.fastq1),
                                readset.fastq2
                            )
                        ],
                        name="fastq_rename.R2.I2." + readset.name + "." + self.run_id + "." + lane,
                        samples=[readset.sample]
                    )
                )
                cleanjob_deps.extend([readset.fastq2, readset.index_fastq2])
            elif ''.join(i for i in self.mask if not i.isdigit()) == "Y,Y,I,Y":
                postprocessing_jobs.append(
                    concat_jobs(
                        [
                            bash.mv(
                                readset.index_fastq1,
                                readset.index_fastq2
                            ),
                            bash.mv(
                                readset.fastq2,
                                readset.index_fastq1
                            ),
                            bash.mv(
                                re.sub("_R1_", "_R3_", readset.fastq1),
                                readset.fastq2
                            )
                        ],
                        name="fastq_rename.R2.I1.I2." + readset.name + "." + self.run_id + "." + lane,
                        samples=[readset.sample]
                    )
                )
                cleanjob_deps.extend([readset.fastq2, readset.index_fastq1, readset.index_fastq2])

        return bcl2fastq_outputs, postprocessing_jobs, cleanjob_deps

    def generate_demuxfastqs_outputs(self, lane):
        demuxfastqs_outputs = []
        postprocessing_jobs = []
        cleanjob_deps = []
        output_dir = self.output_dirs[lane][f"Unaligned.{lane}_directory"]

        index_lengths = self.get_smallest_index_length(lane)
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
                    samples=[readset.sample]
                )
            )
            cleanjob_deps.extend(outputs)

            # Processing "raw R2" fastq (also contains barcode sequences) :
            #   convert headers from MGI to Illumina format
            #   while extracting I1 and I2 to build clean R2, I1 and R2 fastq
            #   using zcat and awk
            if self.is_paired_end[lane]:
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
                        samples=[readset.sample]
                    )
                )
                cleanjob_deps.extend(outputs)

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
        cleanjob_deps.extend(outputs)

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
            if self.is_dual_index[lane]:
                unaligned_i2 = os.path.join(output_dir, "Undetermined_S0_L00" + lane + "_I2_001.fastq.gz")
                outputs.append(unaligned_i2)
                unexpected_barcode_counts_i2 = re.sub(".fastq.gz", ".counts.txt", unaligned_i2)

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
            cleanjob_deps.extend(outputs)

        if unaligned_i1:
            index1length = index_lengths[1] if unaligned_i2 else index_lengths[0]
            if index1length:
                postprocessing_jobs.append(
                    pipe_jobs(
                        [
                            bash.cat(
                                unaligned_i1,
                                None,
                                zip=True
                            ),
                            bash.awk(
                                None,
                                unexpected_barcode_counts_i1,
                                f"'NR%4==2 {{print substr($0,1,{index1length})}}' | sort | uniq -c | sort -nr"
                            )
                        ],
                        name="fastq_countbarcodes.I1.unmatched." + self.run_id + "." + lane,
                        samples=self.samples[lane]
                    )
                )
                cleanjob_deps.append(unexpected_barcode_counts_i1)
        if unaligned_i2:
            index1length = index_lengths[1]
            index2length = index_lengths[0]
            if index2length:
                postprocessing_jobs.append(
                    pipe_jobs(
                        [
                            bash.cat(
                                unaligned_i2,
                                None,
                                zip=True
                            ),
                            bash.awk(
                                None,
                                unexpected_barcode_counts_i2,
                                f"'NR%4==2 {{print substr($0,1,{index2length})}}' | sort | uniq -c | sort -nr"
                            )
                        ],
                        name="fastq_countbarcodes.I2.unmatched." + self.run_id + "." + lane,
                        samples=self.samples[lane]
                    )
                )
                cleanjob_deps.append(unexpected_barcode_counts_i2)
            unexpected_barcode_counts_i1i2 = re.sub("_I1_", "_I1I2_", unexpected_barcode_counts_i1)
            postprocessing_jobs.append(
                concat_jobs(
                    [
                        bash.paste(
                            None,
                            unexpected_barcode_counts_i1i2,
                            options=f"-d '' <(zcat {unaligned_i1} | awk 'NR%4==2 {{print substr($0,1,{index1length})}}') <(zcat {unaligned_i2} | awk 'NR%4==2 {{print substr($0,1,{index2length})}}') | sort | uniq -c | sort -nr"
                        )
                    ],
                    input_dependency=[
                        unaligned_i1,
                        unaligned_i2
                    ],
                    name="fastq_countbarcodes.I1I2.unmatched." + self.run_id + "." + lane,
                    samples=self.samples[lane]
                )
            )
            cleanjob_deps.append(unexpected_barcode_counts_i1)

        #candidate_input_files=[[unexpected_barcode_counts_i1i2], [unexpected_barcode_counts_i1]]
        #unexpected_barcode_counts = self.select_input_files(candidate_input_files)
        #unexpected_barcode_matches = re.sub(".counts.txt", ".match_table.tsv", unexpected_barcode_counts)
        if unaligned_i2:
            unexpected_barcode_counts = unexpected_barcode_counts_i1i2
        else:
            unexpected_barcode_counts = unexpected_barcode_counts_i1
        unexpected_barcode_matches = re.sub(".counts.txt", ".match_table.tsv", unexpected_barcode_counts)

        job = run_processing_tools.match_undetermined_barcodes(
                unexpected_barcode_counts,
                unexpected_barcode_matches
                )
        job.name = 'fastq_match_undetermined_barcodes.' + self.run_id + "." + lane
        samples = self.samples[lane]
        postprocessing_jobs.append(job)

        return demuxfastqs_outputs, postprocessing_jobs, cleanjob_deps

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
                    r1_out=readset.fastq1 if readset else os.path.join(self.output_dirs[lane][f"Unaligned.{lane}_directory"], "Undetermined_S0_L00" + lane + "_R1_001.fastq.gz"),
                    i1_out=readset.index_fastq1 if readset else os.path.join(self.output_dirs[lane][f"Unaligned.{lane}_directory"], "Undetermined_S0_L00" + lane + "_I1_001.fastq.gz"),
                    i2_out=readset.index_fastq2 if readset else os.path.join(self.output_dirs[lane][f"Unaligned.{lane}_directory"], "Undetermined_S0_L00" + lane + "_I2_001.fastq.gz")
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
                    r1_out=readset.fastq1 if readset else os.path.join(self.output_dirs[lane][f"Unaligned.{lane}_directory"], "Undetermined_S0_L00" + lane + "_R1_001.fastq.gz"),
                    i1_out=readset.index_fastq1 if readset else os.path.join(self.output_dirs[lane][f"Unaligned.{lane}_directory"], "Undetermined_S0_L00" + lane + "_I1_001.fastq.gz")
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
                r2_out=readset.fastq2 if readset else os.path.join(self.output_dirs[lane][f"Unaligned.{lane}_directory"], "Undetermined_S0_L00" + lane + "_R2_001.fastq.gz"),
                i1_out=readset.index_fastq1 if readset else os.path.join(self.output_dirs[lane][f"Unaligned.{lane}_directory"], "Undetermined_S0_L00" + lane + "_I1_001.fastq.gz"),
                i2_out=readset.index_fastq2 if readset else os.path.join(self.output_dirs[lane][f"Unaligned.{lane}_directory"], "Undetermined_S0_L00" + lane + "_I2_001.fastq.gz")
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
                r2_out=readset.fastq2 if readset else os.path.join(self.output_dirs[lane][f"Unaligned.{lane}_directory"], "Undetermined_S0_L00" + lane + "_R2_001.fastq.gz"),
                i1_out=readset.index_fastq1 if readset else os.path.join(self.output_dirs[lane][f"Unaligned.{lane}_directory"], "Undetermined_S0_L00" + lane + "_I1_001.fastq.gz")
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
        return [ RunInfoRead(int(r.get("Number")), int(r.get("NumCycles")), r.get("IsIndexedRead") == "Y", r.get("IsReverseComplement") == "Y") for r in reads.iter('Read') ]
    
    def get_run_parameter_file(self, run_dir):
        """
        Find the Illumina run parameter XML file within the run folder : could be 'runParameters.xml' or 'RunParameters.xml'
        """
        param_files = [filename for filename in os.listdir(run_dir) if re.match("RunParameters.xml", filename, flags=re.IGNORECASE)]
        if len(param_files) > 1:
            _raise(SanitycheckError("More than one run parameter XML file found in " + run_dir + ":\n\t" + "\t".join(param_files)))
        elif len(param_files) == 0:
            _raise(SanitycheckError("No run parameter XML file found in " + run_dir + "..."))
        return os.path.join(run_dir, param_files[0])

    def get_sbs_consumable_version(self):
        """
        Parse the Illumina RunParameters.xml file for the 'SbsConsumableVersion'
        """
        if self.args.type == 'illumina':
            return Xml.parse(self.get_run_parameter_file(self.run_dir)).findall('.//SbsConsumableVersion')[0].text

    def load_readsets(self, lane):
        """
        Parse the sample sheet and return a list of readsets.
        """
        seqtype = "hiseqx" if (self.seqtype == "novaseq" and self.sbs_consumable_version == '3') else self.seqtype
        if is_json(self.readset_file):
            return parse_freezeman_readset_file(
                self.readset_file,
                self.run_id,
                "PAIRED_END" if self.is_paired_end[lane] else "SINGLE_END",
                lane,
                seqtype,
                int(self.read1cycles[lane]),
                int(self.read2cycles[lane]),
                int(self.index1cycles[lane]),
                int(self.index2cycles[lane]),
                self.output_dirs,
                self.args.type
            )
        else:
            return parse_clarity_readset_file(
                self.readset_file,
                self.run_id,
                "PAIRED_END" if self.is_paired_end[lane] else "SINGLE_END",
                lane,
                seqtype,
                int(self.read1cycles[lane]),
                int(self.read2cycles[lane]),
                int(self.index1cycles[lane]),
                int(self.index2cycles[lane]),
                self.output_dirs,
                self.args.type
            )

    def submit_jobs(self):
        super(RunProcessing, self).submit_jobs()

    @property
    def steps(self):
        # [ illumina, mgig400, mgit7 ]
        return [
            [
                self.index,
                self.fastq,
                self.fastp,
                self.blast,
                self.align,
                self.picard_mark_duplicates,
                self.check_sample_mixup,
                self.metrics,
                self.md5,
                self.report,
                self.copy,
                self.final_notification
            ],
            [
                self.fastq,
                self.fastp,
                self.blast,
                self.align,
                self.picard_mark_duplicates,
                self.check_sample_mixup,
                self.metrics,
                self.md5,
                self.report,
                self.copy,
                self.final_notification
            ],
            [
                self.basecall,
                self.fastq,
                self.fastp,
                self.blast,
                self.align,
                self.picard_mark_duplicates,
                self.check_sample_mixup,
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

def is_json(filepath):
    """
    Checks whether a file is a JSON file or not.
    Returns True or False
    """
    with open(filepath) as f:
        if f.read(1) in '{[':
            return True
        else:
            return False

if __name__ == '__main__':

    argv = sys.argv
    if '--wrap' in argv:
        utils.container_wrapper_argparse(argv)
    else:
        RunProcessing(protocol=['illumina', 'mgig400', 'mgit7'])

