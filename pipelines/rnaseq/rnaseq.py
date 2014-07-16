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

from bio import bedtools
from bio import metrics
from bio import picard
from bio import samtools
from bio import tophat
from pipelines.illumina import illumina

log = logging.getLogger(__name__)

class RnaSeq(illumina.Illumina):

    def trim_metrics(self):
        # Transform pipeline 'PAIRED_END' or 'SINGLE_END' run_type to 'paired' or 'single' parameter
        run_type = re.sub("_END$", "", self.run_type).lower()

        job = metrics.merge_trimmomatic_stats("trim.stats.csv", "trim", os.path.join("metrics", "trimming.stats"), run_type)
        job.input_files = [os.path.join("trim", readset.sample.name, readset.name + ".trim.stats.csv") for readset in self.readsets]
        job.name = "trim_metrics"
        return [job]

    def tophat(self):
        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_directory = os.path.join("alignment", readset.sample.name)
            readset_alignment_directory = os.path.join(alignment_directory, readset.name)

            if readset.run_type == "PAIRED_END":
                fastq1 = trim_file_prefix + "pair1.fastq.gz"
                fastq2 = trim_file_prefix + "pair2.fastq.gz"
            elif readset.run_type == "SINGLE_END":
                fastq1 = trim_file_prefix + "single.fastq.gz"
                fastq2 = None
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            job = tophat.tophat(
                fastq1,
                fastq2,
                readset_alignment_directory,
                rg_id=readset.name,
                rg_sample=readset.sample.name,
                rg_library=readset.library,
                rg_platform_unit=readset.run + "_" + readset.lane,
                rg_platform=config.param('tophat', 'platform'),
                rg_center=config.param('tophat', 'TBInstitution'),
            )

            # If this readset is unique for this sample, further BAM merging is not necessary.
            # Thus, create a sample BAM symlink to the readset BAM.
            if len(readset.sample.readsets) == 1:
                sample_bam = os.path.join(alignment_directory, readset.sample.name + ".sorted.bam")
                job.command += " && \\\nln -s " + os.path.join(readset_alignment_directory, "accepted_hits.bam") + " " + sample_bam
                job.output_files.append(sample_bam)

            job.name = "tophat." + readset.name
            jobs.append(job)

        return jobs

    def picard_merge_sam_files(self):
        jobs = []
        for sample in self.samples:
            # Skip samples with one readset only, since symlink has been created at align step
            if len(sample.readsets) > 1:
                alignment_directory = os.path.join("alignment", sample.name)
                inputs = [os.path.join(alignment_directory, readset.name + ".sorted.bam") for readset in sample.readsets]
                output = os.path.join(alignment_directory, sample.name + ".sorted.bam")

                job = picard.merge_sam_files(inputs, output)
                job.name = "picard_merge_sam_files." + sample.name
                jobs.append(job)
        return jobs

    def picard_reorder_sam(self):
        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name)

            job = picard.reorder_sam(
                alignment_file_prefix + ".sorted.bam",
                alignment_file_prefix + ".merged.karyotypic.bam"
            )
            job.name = "picard_reorder_sam." + sample.name
            jobs.append(job)
        return jobs

    def picard_mark_duplicates(self):
        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".merged.")

            job = picard.mark_duplicates(
                [alignment_file_prefix + "karyotypic.bam"],
                alignment_file_prefix + "mdup.bam",
                alignment_file_prefix + "mdup.metrics"
            )
            job.name = "picard_mark_duplicates." + sample.name
            jobs.append(job)
        return jobs

    def rnaseqc(self):

        sample_file = os.path.join("alignment", "rnaseqc.samples.txt")
        output_directory = os.path.join("metrics", "rnaseqRep")

        job = metrics.rnaseqc(sample_file, output_directory, self.run_type == "SINGLE_END")

        sample_file_job = Job([os.path.join("alignment", sample.name, sample.name + ".merged.mdup.bam") for sample in self.samples], [sample_file])
        project_name = config.param('DEFAULT', 'projectName')
        job.command = """\
mkdir -p {output_directory} && \\
echo \\"Sample\tBamFile\tNote
{input_bams}\\" \
> {sample_file} && \\
{job.command} && \\
zip -r {output_directory}.zip {output_directory}""".format(
            output_directory=output_directory,
            input_bams=" \\\n".join(["\t".join([sample.name, os.path.join("alignment", sample.name, sample.name + ".merged.mdup.bam"), project_name]) for sample in self.samples]),
            sample_file=sample_file,
            job=job
        )

        job.input_files.extend([os.path.join("alignment", sample.name, sample.name + ".merged.mdup.bam") for sample in self.samples])
        job.output_files.extend([sample_file, output_directory + ".zip"])
        job.name = "rnaseqc"
        return [job]

    def wiggle(self):
        jobs = []

        for sample in self.samples:
            bam_file_prefix = os.path.join("alignment", sample.name, sample.name + ".merged.mdup.")
            input_bam = bam_file_prefix + "bam"
            bed_graph_prefix = os.path.join("tracks", sample.name, sample.name)
            big_wig_prefix = os.path.join("tracks", "bigWig", sample.name)

            if config.param('tophat', 'strandInfo') != 'fr-unstranded':
                input_bam_f1 = bam_file_prefix + "tmp1.forward.bam"
                input_bam_f2 = bam_file_prefix + "tmp2.forward.bam"
                input_bam_r1 = bam_file_prefix + "tmp1.reverse.bam"
                input_bam_r2 = bam_file_prefix + "tmp2.reverse.bam"
                output_bam_f = bam_file_prefix + "forward.bam"
                output_bam_r = bam_file_prefix + "reverse.bam"

                bam_f_job = concat_jobs([
                    samtools.view(input_bam, input_bam_f1, "-bh -F 256 -f 81"),
                    samtools.view(input_bam, input_bam_f2, "-bh -F 256 -f 161"),
                    picard.merge_sam_files([input_bam_f1, input_bam_f2], output_bam_f),
                    Job(command="rm " + input_bam_f1 + " " + input_bam_f2)
                ])
                bam_f_job.name = "wiggle." + sample.name + ".forward_strandspec"

                bam_r_job = concat_jobs([
                    Job(command="mkdir -p " + os.path.join("tracks", sample.name) + " " + os.path.join("tracks", "bigWig")),
                    samtools.view(input_bam, input_bam_r1, "-bh -F 256 -f 97"),
                    samtools.view(input_bam, input_bam_r2, "-bh -F 256 -f 145"),
                    picard.merge_sam_files([input_bam_r1, input_bam_r2], output_bam_r),
                    Job(command="rm " + input_bam_r1 + " " + input_bam_r2)
                ])
                bam_r_job.name = "wiggle." + sample.name + ".reverse_strandspec"

                jobs.extend([bam_f_job, bam_r_job])

                outputs = [
                    [bed_graph_prefix + ".forward.bedGraph", big_wig_prefix + ".forward.bw"],
                    [bed_graph_prefix + ".reverse.bedGraph", big_wig_prefix + ".reverse.bw"],
                ]
            else:
                outputs = [[bed_graph_prefix + ".bedGraph", big_wig_prefix + ".bw"]]

            for bed_graph_output, big_wig_output in outputs:
                job = concat_jobs([
                    Job(command="mkdir -p " + os.path.join("tracks", sample.name) + " " + os.path.join("tracks", "bigWig")),
                    bedtools.graph(input_bam, bed_graph_output, big_wig_output)
                ])
                job.name = "wiggle." + re.sub(".bedGraph", "", os.path.basename(bed_graph_output))
                jobs.append(job)

        return jobs

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.trim_metrics,
            self.tophat,
            self.picard_merge_sam_files,
            self.picard_reorder_sam,
            self.picard_mark_duplicates,
            self.rnaseqc,
            self.wiggle
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
