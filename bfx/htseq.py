#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def htseq_count(input, gtf, output, options="", stranded="no"):
    return Job(
        [input],
        [output],
        [['htseq_count', 'module_python']],
        command="""\
htseq-count {options} \\
  --stranded={stranded} \\
  --format=bam \\
  {input} \\
  {gtf}{output}""".format(
        options=options,
        stranded=stranded,
        input=input,
        gtf=gtf,
        output=" \\\n  > " + output if output else ""
        )
    )
