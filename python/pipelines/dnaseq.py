#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections
import logging
import os
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(sys.argv[0])))

# MUGQIC Modules
from core.pipeline import *
from bio.readset import *
from bio.trimmomatic import *
from bio.bwa import *

log = logging.getLogger(__name__)

class DnaSeq(Pipeline):

    @property
    def readsets(self):
        return self._readsets

    @property
    def samples(self):
        return self._samples

    def trim(self, readset):
        trim_file_prefix = "trim/" + readset.sample.name + "/" + readset.name + ".trim."
        return trimmomatic(
            self.config,
            readset.fastq1,
            readset.fastq2,
            trim_file_prefix + "pair1.fastq.gz",
            trim_file_prefix + "single1.fastq.gz",
            trim_file_prefix + "pair2.fastq.gz",
            trim_file_prefix + "single2.fastq.gz",
            None,
            readset.quality_offset,
            trim_file_prefix + "out",
            trim_file_prefix + "stats.csv"
        )

    def mem(self, readset):
        trim_file_prefix = "trim/" + readset.sample.name + "/" + readset.name + ".trim."

        rg_platform_unit = readset.run + "_" + readset.lane
        rg_id = readset.library + "_" + rg_platform_unit

        read_group = "'@RG\tID:" + rg_id + "\tSM:" + readset.sample.name + "\tLB:" + readset.library + "\tPU:run" + rg_platform_unit + "\tCN:" + self.config.param('mem', 'bwaInstitution') + "\tPL:Illumina'"

        return mem(
            self.config,
            trim_file_prefix + "pair1.fastq.gz",
            trim_file_prefix + "pair2.fastq.gz",
            "align/" + readset.sample.name + "/" + readset.name + ".sam",
            read_group
        )

    @property
    def step_dict_map(self):
        return [
            {"name": self.trim, "loop": self.readsets},
            {"name": self.mem, "loop": self.readsets}
        ]

    def __init__(self):
        # Initialize attributes to avoid AttributeError in default_argparser(self.step_dict_map)
        self._readsets = []
        self._samples = []

        argparser = PipelineArgumentParser(self.step_dict_map)
        # Add pipeline specific arguments
        argparser.add_argument("-r", "--readsets", help="readset file", type=file, required=True)
        args = argparser.parse_args()

        self._readsets = parse_readset_file(args.readsets.name)
        # Retrieve unique samples from their readsets, removing duplicates
        self._samples = list(collections.OrderedDict.fromkeys([readset.sample for readset in self._readsets]))

        Pipeline.__init__(self, args)
        
DnaSeq().submit_jobs()
