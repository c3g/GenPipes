#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def htseq_count(input, output, options="", stranded="no"):
    job = Job([input], [output], [['htseq_count', 'moduleVersion.python']])

    job.command = \
"""htseq-count {options} \\
  --stranded={stranded} \\
  {input}{output}""".format(
        options=options,
        stranded=stranded,
        input=input,
        output=" \\\n  > " + output if output else ""
    )

    return job
