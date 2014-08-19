#!/usr/bin/env python

# Python Standard Modules
import logging
import os
import re
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.job import *
from core.pipeline import *
from bio.design import *
from bio.readset import *

from bio import metrics
from bio import picard
from bio import trimmomatic

log = logging.getLogger(__name__)

# Abstract pipeline gathering common features of all Illumina sequencing pipelines (readsets, trimming, etc.)
# Specific steps must be defined in Illumina children pipelines.
class Illumina(Pipeline):

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = parse_readset_file(self.args.readsets.name)
        return self._readsets

    @property
    def samples(self):
        if not hasattr(self, "_samples"):
            self._samples = list(collections.OrderedDict.fromkeys([readset.sample for readset in self.readsets]))
        return self._samples

    @property
    def run_type(self):
        run_types = [readset.run_type for readset in self.readsets]
        if len(set(run_types)) == 1 and re.search("^(PAIRED|SINGLE)_END$", run_types[0]):
            return run_types[0]
        else:
            raise Exception("Error: readset run types " + ",".join(["\"" + run_type + "\"" for run_type in run_types]) +
            " are invalid (should be all PAIRED_END or all SINGLE_END)!")

    @property
    def contrasts(self):
        if not hasattr(self, "_contrasts"):
            if self.args.design:
                self._contrasts = parse_old_design_file(self.args.design.name, self.samples)
            else:
                raise Exception("Error: missing '--design' option!")
        return self._contrasts

    def picard_sam_to_fastq(self):
        jobs = []
        for readset in self.readsets:
            if readset.bam and not readset.fastq1:
                if readset.run_type == "PAIRED_END":
                    readset.fastq1 = re.sub("\.bam$", ".pair1.fastq.gz", readset.bam)
                    readset.fastq2 = re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)
                elif readset.run_type == "SINGLE_END":
                    readset.fastq1 = re.sub("\.bam$", ".single.fastq.gz", readset.bam)
                else:
                    raise Exception("Error: run type \"" + readset.run_type +
                    "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

                job = picard.sam_to_fastq(readset.bam, readset.fastq1, readset.fastq2)
                job.name = "picard_sam_to_fastq." + readset.name
                jobs.append(job)
        return jobs

    def trimmomatic(self):
        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            trim_log = trim_file_prefix + "log"
            trim_stats = trim_file_prefix + "stats.csv"
            if readset.run_type == "PAIRED_END":
                if readset.bam and not readset.fastq1:
                    readset.fastq1 = re.sub("\.bam$", ".pair1.fastq.gz", readset.bam)
                    readset.fastq2 = re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)
                job = trimmomatic.trimmomatic(
                    readset.fastq1,
                    readset.fastq2,
                    trim_file_prefix + "pair1.fastq.gz",
                    trim_file_prefix + "single1.fastq.gz",
                    trim_file_prefix + "pair2.fastq.gz",
                    trim_file_prefix + "single2.fastq.gz",
                    None,
                    readset.quality_offset,
                    trim_log
                )
            elif readset.run_type == "SINGLE_END":
                if readset.bam and not readset.fastq1:
                    readset.fastq1 = re.sub("\.bam$", ".single.fastq.gz", readset.bam)
                job = trimmomatic.trimmomatic(
                    readset.fastq1,
                    None,
                    None,
                    None,
                    None,
                    None,
                    trim_file_prefix + "single.fastq.gz",
                    readset.quality_offset,
                    trim_log
                )
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            job.name = "trimmomatic." + readset.name
            jobs.append(job)
        return jobs

    def merge_trimmomatic_stats(self):
        merge_trim_stats = os.path.join("metrics", "trimming.stats")
        job = Job(command="rm -f " + merge_trim_stats)
        for readset in self.readsets:
            trim_log = os.path.join("trim", readset.sample.name, readset.name + ".trim.log")
            if readset.run_type == "PAIRED_END":
                perl_command = "perl -pe 's/^Input Read Pairs: (\\\\d+).*Both Surviving: (\\\\d+).*Forward Only Surviving: (\\\\d+).*$/{readset.sample.name}\t{readset.name}\t\\\\1\t\\\\2\t\\\\3/'".format(readset=readset)
            elif readset.run_type == "SINGLE_END":
                perl_command = "perl -pe 's/^Input Reads: (\\\\d+).*Surviving: (\\\\d+).*$/{readset.sample.name}\t{readset.name}\t\\\\1\t\\\\2\t\\\\2/'".format(readset=readset)

            job = concat_jobs([
                job,
                Job(
                    [trim_log],
                    [merge_trim_stats],
                    command="""\
grep ^Input {trim_log} | \\
{perl_command} \\
  >> {merge_trim_stats}""".format(
                        trim_log=trim_log,
                        perl_command=perl_command,
                        merge_trim_stats=merge_trim_stats
                    )
                )
            ], name="merge_trimmomatic_stats")

        return [job]

    def __init__(self):
        # Add pipeline specific arguments
        self.argparser.add_argument("-r", "--readsets", help="readset file", type=file, required=True)

        super(Illumina, self).__init__()
