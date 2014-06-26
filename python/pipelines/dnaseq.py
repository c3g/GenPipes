#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(sys.argv[0])))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bio.readset import *
from bio.trimmomatic import *
from bio.bwa import *
from bio.gatk import *
from bio.picard import *
from bio.sequence_dictionary import *

log = logging.getLogger(__name__)

class DnaSeq(Pipeline):

    @property
    def readsets(self):
        return self._readsets

    @property
    def samples(self):
        return self._samples

    @property
    def sequence_dictionary(self):
        return self._sequence_dictionary

    def sam_to_fastq(self, readset):
        if readset.bam and not readset.fastq1:
            if readset.run_type == "PAIRED_END":
                readset.fastq1 = re.sub("\.bam$", ".pair1.fastq.gz", readset.bam)
                readset.fastq2 = re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)
            elif readset.run_type == "SINGLE_END":
                fastq1 = re.sub("\.bam$", ".single.fastq.gz", readset.bam)
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            return sam_to_fastq(readset.bam, readset.fastq1, readset.fastq2)

    def trim(self, readset):
        trim_file_prefix = "trim/" + readset.sample.name + "/" + readset.name + ".trim."
        return trimmomatic(
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

    def bwa_mem_sort_sam(self, readset):
        trim_file_prefix = "trim/" + readset.sample.name + "/" + readset.name + ".trim."

        rg_platform_unit = readset.run + "_" + readset.lane
        rg_id = readset.library + "_" + rg_platform_unit

        read_group = "'@RG\tID:" + rg_id + "\tSM:" + readset.sample.name + "\tLB:" + readset.library + "\tPU:run" + rg_platform_unit + "\tCN:" + config.param('mem', 'bwaInstitution') + "\tPL:Illumina'"

        bwa_job = mem(
            trim_file_prefix + "pair1.fastq.gz",
            trim_file_prefix + "pair2.fastq.gz",
            None,
            read_group
        )

        sort_sam_job = sort_sam(
            "/dev/stdin",
            "align/" + readset.sample.name + "/" + readset.name + ".sorted.bam",
            "coordinate"
        )

        return pipe_jobs([bwa_job, sort_sam_job])

    def merge_readsets(self, sample):
        align_file_prefix = "align/" + sample.name + "/"
        inputs = [align_file_prefix + readset.name + ".sorted.bam" for readset in sample.readsets]
        output = align_file_prefix + sample.name + ".sorted.bam"

        return merge_sam_files(inputs, output)

    def mark_duplicates(self, sample):
        align_file_prefix = "align/" + sample.name + "/" + sample.name + ".sorted."
        input = align_file_prefix + "bam"
        output = align_file_prefix + "dup.bam"
        metrics_file = align_file_prefix + "dup.metrics"

        return mark_duplicates([input], output, metrics_file)

    def collect_multiple_metrics(self, sample):
        align_file_prefix = "align/" + sample.name + "/" + sample.name + ".sorted."
        input = align_file_prefix + "bam"
        output = align_file_prefix + "dup.metrics"

        return collect_multiple_metrics(input, output)


    def recalibration(self, sample):
        align_file_prefix = "align/" + sample.name + "/" + sample.name + ".sorted.dup."
        input = align_file_prefix + "bam"
        print_reads_output = align_file_prefix + "recal.bam"
        base_recalibrator_output = align_file_prefix + "recalibration_report.grp"

        return concat_jobs([
            base_recalibrator(input, base_recalibrator_output),
            print_reads(input, print_reads_output, base_recalibrator_output)
        ])

    @property
    def step_dict_map(self):
        return [
            {"name": self.sam_to_fastq, "loop": self.readsets},
            {"name": self.trim, "loop": self.readsets},
            {"name": self.bwa_mem_sort_sam, "loop": self.readsets},
            {"name": self.merge_readsets, "loop": self.samples},
            {"name": self.mark_duplicates, "loop": self.samples},
            {"name": self.collect_multiple_metrics, "loop": self.samples},
            {"name": self.recalibration, "loop": self.samples}
        ]

    def __init__(self):
        # Initialize attributes to avoid AttributeError in default_argparser(self.step_dict_map)
        self._readsets = []
        self._samples = []

        argparser = PipelineArgumentParser(self.step_dict_map)
        # Add pipeline specific arguments
        argparser.add_argument("-r", "--readsets", help="readset file", type=file, required=True)
        args = argparser.parse_args()

        # Create sequence dictionary
        self._sequence_dictionary = parse_sequence_dictionary_file(config.param('DEFAULT', 'referenceSequenceDictionary', type='filepath'))

        # Create readsets
        self._readsets = parse_readset_file(args.readsets.name)

        # Retrieve unique samples from their readsets, removing duplicates
        self._samples = list(collections.OrderedDict.fromkeys([readset.sample for readset in self._readsets]))

        Pipeline.__init__(self, args)
        
DnaSeq().submit_jobs()
