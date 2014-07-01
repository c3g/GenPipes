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

def sort(input_bam, output_prefix):
    job = Job([input_bam], [output_prefix + ".bam"], [['samtools_sort', 'moduleVersion.samtools']])

    job.command = \
"""samtools sort {extra_sort_flags} {input_bam} {output_prefix}""".format(
        extra_sort_flags=config.param('samtools_sort', 'extraSortFlags', required=False),
        input_bam=input_bam,
        output_prefix=output_prefix
    )

    return job
