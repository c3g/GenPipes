#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def index(input, output):
    job = Job([input], [output], [['samtools_index', 'module_samtools']])

    job.command = \
"""samtools index \\
  {input}""".format(
        input=input
    )

    return job

def faidx(input):
    job = Job([input], [input + ".fai"], [['samtools_index', 'module_samtools']])

    job.command = \
"""samtools faidx \\
  {input}""".format(
        input=input
    )

    return job

def flagstat(input, output):
    job = Job([input], [output], [['samtools_flagstat', 'module_samtools']])

    job.command = \
"""samtools flagstat \\
  {input} \\
  > {output}""".format(
        input=input,
        output=output
    )

    return job

def mpileup(input_bams, output, other_options="", region=None, pair_calling=False):
    job = Job(input_bams, [output], [['samtools_mpileup', 'module_samtools']])

    job.command = \
"""samtools mpileup {other_options} \\
  -f {reference_fasta}{region}{input_bams}{output}""".format(
        other_options=other_options,
        reference_fasta=config.param('samtools_mpileup', 'genome_fasta', type='filepath'),
        region=" \\\n  -r " + region if region else "",
        input_bams="".join([" \\\n  " + input_bam for input_bam in input_bams]),
        output=" \\\n  > " + output if output else ""
    )

    return job

def sort(input_bam, output_prefix, sort_by_name=False):
    job = Job([input_bam], [output_prefix + ".bam"], [['samtools_sort', 'module_samtools']])

    other_options=" -n " if sort_by_name else ""
    
    job.command = \
"""samtools sort {other_options} {input_bam} {output_prefix}""".format(
        other_options=other_options + config.param('samtools_sort', 'other_options', required=False),
        input_bam=input_bam,
        output_prefix=output_prefix
    )

    return job

def view(input, output=None, options=""):
    job = Job([input], [output], [['samtools_view', 'module_samtools']])

    job.command = \
"""samtools view {options} \\
  {input}{output}""".format(
        options=options,
        input=input,
        output=" \\\n  > " + output if output else ""
    )

    return job

def bcftools_cat(inputs, output):
    job = Job(inputs, [output], [['bcftools_cat', 'module_samtools']])

    job.command = \
"""bcftools cat \\
  {inputs}{output}""".format(
        inputs=" \\\n  ".join(inputs),
        output=" \\\n  > " + output if output else ""
    )

    return job

def bcftools_view(input, output, options="", pair_calling=False):
    job = Job([input], [output], [['bcftools_view', 'module_samtools']])

    job.command = \
"""bcftools view {pair_calling} {options} \\
  {input}{output}""".format(
        options=options,
        pair_calling="-T pair" if pair_calling else "",
        input=input,
        output=" \\\n  > " + output if output else ""
    )

    return job
