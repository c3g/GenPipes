#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.job import *

def mem(
    config,
    in1fastq,
    in2fastq=None,
    out_sam=None,
    read_group=None
    ):

    job = Job(config, [in1fastq, in2fastq], [out_sam], [["mem", "moduleVersion.bwa"]])

    job.command = ""
    if out_sam:
        job.command += "mkdir -p " + os.path.dirname(out_sam) + " && \\\n"

    idxbase = config.param('mem', 'bwaRefIndex', type='filepath')
    extra_flags = config.param('mem', 'bwaExtraFlags')

    job.command += \
"""bwa mem {extra_flags}{read_group} \\
  {idxbase} \\
  {in1fastq}{in2fastq}{out_sam}""".format(
        extra_flags=" \\\n  " + extra_flags if extra_flags else "",
        read_group=" \\\n  -R " + read_group if read_group else "",
        idxbase=idxbase,
        in1fastq=in1fastq,
        in2fastq=" \\\n  " + in2fastq if in2fastq else "",
        out_sam=" \\\n  > " + out_sam if out_sam else ""
    )

    return job
