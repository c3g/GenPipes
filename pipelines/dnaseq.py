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

                job = sam_to_fastq(readset.bam, readset.fastq1, readset.fastq2)
                job.name = "sam_to_fastq." + readset.name
                jobs.append(job)
        return jobs

    def trim(self):
        jobs = []
        for readset in self.readsets:
            trim_file_prefix = "trim/" + readset.sample.name + "/" + readset.name + ".trim."
            job = trimmomatic(
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

    def bwa_mem_sort_sam(self):
        jobs = []
        for readset in self.readsets:
            trim_file_prefix = "trim/" + readset.sample.name + "/" + readset.name + ".trim."
            alignment_directory = "alignment/" + readset.sample.name + "/"

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
                alignment_directory + readset.name + ".sorted.bam",
                "coordinate"
            )

            job = pipe_jobs([bwa_job, sort_sam_job])

            # Create alignment directory (not done by default by BWA mem or Picard SortSam)
            job.command = "mkdir -p " + alignment_directory + " && \\\n" + job.command

            job.name = "bwa_mem_sort_sam." + readset.name
            jobs.append(job)
        return jobs

    def merge_readsets(self):
        jobs = []
        for sample in self.samples:
            align_file_prefix = "alignment/" + sample.name + "/"
            inputs = [align_file_prefix + readset.name + ".sorted.bam" for readset in sample.readsets]
            output = align_file_prefix + sample.name + ".sorted.bam"

            job = merge_sam_files(inputs, output)
            job.name = "merge_readsets." + sample.name
            jobs.append(job)
        return jobs

    def mark_duplicates(self):
        jobs = []
        for sample in self.samples:
            align_file_prefix = "alignment/" + sample.name + "/" + sample.name + ".sorted."
            input = align_file_prefix + "bam"
            output = align_file_prefix + "dup.bam"
            metrics_file = align_file_prefix + "dup.metrics"

            job = mark_duplicates([input], output, metrics_file)
            job.name = "mark_duplicates." + sample.name
            jobs.append(job)
        return jobs

    def collect_multiple_metrics(self):
        jobs = []
        for sample in self.samples:
            align_file_prefix = "alignment/" + sample.name + "/" + sample.name + ".sorted."
            input = align_file_prefix + "bam"
            output = align_file_prefix + "dup.metrics"

            job = collect_multiple_metrics(input, output)
            job.name = "collect_multiple_metrics." + sample.name
            jobs.append(job)
        return jobs

    def recalibration(self):
        jobs = []
        for sample in self.samples:
            align_file_prefix = "alignment/" + sample.name + "/" + sample.name + ".sorted.dup."
            input = align_file_prefix + "bam"
            print_reads_output = align_file_prefix + "recal.bam"
            base_recalibrator_output = align_file_prefix + "recalibration_report.grp"

            job = concat_jobs([
                base_recalibrator(input, base_recalibrator_output),
                print_reads(input, print_reads_output, base_recalibrator_output)
            ])
            job.name = "recalibration." + sample.name
            jobs.append(job)
        return jobs

    def indel_realigner(self):
        jobs = []

        # Create sequence dictionary
        self._sequence_dictionary = parse_sequence_dictionary_file(config.param('DEFAULT', 'referenceSequenceDictionary', type='filepath'))

        nb_realign_jobs = config.param('indel_realigner', 'nbRealignJobs', type='posint')
        if nb_realign_jobs > 50:
            log.warning("Number of realign jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for sample in self.samples:
            align_directory = "alignment/" + sample.name + "/"
            realign_directory = align_directory + "realign/"
            input = align_directory + sample.name + ".sorted.bam"
            if nb_realign_jobs == 1:
                output_prefix = realign_directory + "all"
                job = concat_jobs([
                    realigner_target_creator(input, output_prefix + ".intervals"),
                    indel_realigner(input, output_prefix + ".bam", target_intervals=output_prefix + ".intervals")
                ])
                # Create output directory since it is done by default by GATK tools
                job.command = "mkdir -p " + realign_directory + " && \\\n" + job.command
                job.name = "indel_realigner." + sample.name
                jobs.append(job)
            else:
                unique_sequences_per_job = [sequence['name'] for sequence in self._sequence_dictionary[0:min(nb_realign_jobs, len(self._sequence_dictionary))]]
                for sequence in unique_sequences_per_job:
                    output_prefix = realign_directory + sequence
                    job = concat_jobs([
                        realigner_target_creator(input, output_prefix + ".intervals", intervals=[sequence]),
                        indel_realigner(input, output_prefix + ".bam", target_intervals=output_prefix + ".intervals", intervals=[sequence])
                    ])
                    job.name = "indel_realigner." + sample.name + "." + sequence
                    # Create output directory since it is done by default by GATK tools
                    job.command = "mkdir -p " + realign_directory + " && \\\n" + job.command
                    jobs.append(job)
                output_prefix = realign_directory + "others"
                job = concat_jobs([
                    realigner_target_creator(input, output_prefix + ".intervals", exclude_intervals=unique_sequences_per_job),
                    indel_realigner(input, output_prefix + ".bam", target_intervals=output_prefix + ".intervals", exclude_intervals=unique_sequences_per_job)
                ])
                # Create output directory since it is done by default by GATK tools
                job.command = "mkdir -p " + realign_directory + " && \\\n" + job.command
                job.name = "indel_realigner." + sample.name + ".others"
                jobs.append(job)

        return jobs

    @property
    def steps(self):
        return [
            self.sam_to_fastq,
            self.trim,
            self.bwa_mem_sort_sam,
            self.merge_readsets,
            self.mark_duplicates,
            self.collect_multiple_metrics,
            self.recalibration,
            self.indel_realigner
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
        
DnaSeq().submit_jobs()
