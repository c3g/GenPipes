#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def flagstat(input, output):
    job = Job([input], [output], [['samtools_flagstat', 'moduleVersion.samtools']])

    job.command = \
"""samtools flagstat \\
  {input} \\
  > {output}""".format(
        input=input,
        output=output
    )

    return job

def mpileup(input_bams, output, extra_options="", region=None, pair_calling=False):
    job = Job(input_bams, [output], [['samtools_mpileup', 'moduleVersion.samtools']])

    job.command = \
"""samtools mpileup {extra_options} \\
  -f {reference_fasta}{region}{input_bams}{output}""".format(
        extra_options=extra_options,
        reference_fasta=config.param('samtools_mpileup', 'referenceFasta', type='filepath'),
        region=" \\\n  -r " + region if region else "",
        input_bams="".join([" \\\n  " + input_bam for input_bam in input_bams]),
        output=" \\\n  > " + output if output else ""
    )

    return job

def sort(input_bam, output_prefix):
    job = Job([input_bam], [output_prefix + ".bam"], [['samtools_sort', 'moduleVersion.samtools']])

    job.command = \
"""samtools sort {extra_sort_flags} {input_bam} {output_prefix}""".format(
        extra_sort_flags=config.param('samtools_sort', 'extraSortFlags', required=False),
        input_bam=input_bam,
        output_prefix=output_prefix
    )

    return job

def samtools_view(input, output, options=""):
    job = Job([input], [output], [['samtools_view', 'moduleVersion.samtools']])

    job.command = \
"""samtools view {options} \\
  {input}{output}""".format(
        options=options,
        input=input,
        output=" \\\n  > " + output if output else ""
    )

    return job

def bcftools_cat(inputs, output):
    job = Job(inputs, [output], [['bcftools_cat', 'moduleVersion.samtools']])

    job.command = \
"""bcftools cat \\
  {inputs}{output}""".format(
        inputs=" \\\n  ".join(inputs),
        output=" \\\n  > " + output if output else ""
    )

    return job

def bcftools_view(input, output, options="", pair_calling=False):
    job = Job([input], [output], [['bcftools_view', 'moduleVersion.samtools']])

    job.command = \
"""bcftools view {pair_calling} {options} \\
  {input}{output}""".format(
        options=options,
        pair_calling="-T pair" if pair_calling else "",
        input=input,
        output=" \\\n  > " + output if output else ""
    )

    return job
