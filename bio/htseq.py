#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def htseq_count(input, gtf, output, options="", stranded="no"):
    job = Job([input], [output], [['htseq_count', 'moduleVersion.python']])

    job.command = \
"""htseq-count {options} \\
  --stranded={stranded} \\
  {input}
  {gtf}{output}""".format(
        options=options,
        stranded=stranded,
        input=input,
        gtf=gtf,
        output=" \\\n  > " + output if output else ""
    )

    return job
