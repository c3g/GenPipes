#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(sys.argv[0]))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bio.readset import *

from bio import metrics
from bio import picard
from bio import tophat
from bio import trimmomatic

log = logging.getLogger(__name__)

class RnaSeq(Pipeline):

    @property
    def readsets(self):
        return self._readsets

    @property
    def samples(self):
        return self._samples

    def sam_to_fastq(self):
        jobs = []
        for readset in self.readsets:
            if readset.bam and not readset.fastq1:
                if readset.run_type == "PAIRED_END":
                    readset.fastq1 = re.sub("\.bam$", ".pair1.fastq.gz", readset.bam)
                    readset.fastq2 = re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)
                elif readset.run_type == "SINGLE_END":
                    fastq1 = re.sub("\.bam$", ".single.fastq.gz", readset.bam)
                else:
                    raise Exception("Error: run type \"" + readset.run_type +
                    "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

                job = picard.sam_to_fastq(readset.bam, readset.fastq1, readset.fastq2)
                job.name = "sam_to_fastq." + readset.name
                jobs.append(job)
        return jobs

    def trim(self):
        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            job = trimmomatic.trimmomatic(
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
            job.name = "trim." + readset.name
            jobs.append(job)
        return jobs

    def trim_metrics(self):
        run_types = [readset.run_type for readset in self.readsets]
        if len(set(run_types)) == 1 and re.search("^(PAIRED|SINGLE)_END$", run_types[0]):
            # run_type parameter must be either 'paired' or 'single'
            run_type = re.sub("_END$", "", run_types[0]).lower()
        else:
            raise Exception("Error: readset run types " + ",".join(["\"" + run_type + "\"" for run_type in run_types]) +
            " are invalid (should be all PAIRED_END or all SINGLE_END)!")

        job = metrics.merge_trimmomatic_stats("trim.stats.csv", "trim", os.path.join("metrics", "trimming.stats"), run_type)
        job.input_files = [os.path.join("trim", readset.sample.name, readset.name + ".trim.stats.csv") for readset in self.readsets]
        job.name = "trim_metrics"
        return [job]

    def tophat(self):
        jobs = []
        for readset in self.readsets:

            job = tophat.tophat(
                readset.fastq1,
                readset.fastq2,
                os.path.join("alignment", readset.sample.name, readset.name),
                rg_id=readset.name,
                rg_sample=readset.sample.name,
                rg_library=readset.library,
                rg_platform_unit=readset.run + "_" + readset.lane,
                rg_platform=config.param('tophat', 'platform'),
                rg_center=config.param('tophat', 'TBInstitution'),
            )
            job.name = "tophat." + readset.name
            jobs.append(job)

        return jobs

    @property
    def steps(self):
        return [
            self.sam_to_fastq,
            self.trim,
            self.trim_metrics,
            self.tophat
        ]

    def __init__(self):
        argparser = PipelineArgumentParser(self.steps)
        # Add pipeline specific arguments
        argparser.add_argument("-r", "--readsets", help="readset file", type=file, required=True)
        args = argparser.parse_args()

        # Create readsets
        self._readsets = parse_readset_file(args.readsets.name)

        # Retrieve unique samples from their readsets, removing duplicates
        self._samples = list(collections.OrderedDict.fromkeys([readset.sample for readset in self._readsets]))

        Pipeline.__init__(self, args)
        
RnaSeq().submit_jobs()
