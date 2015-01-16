#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def htseq_count(input, gtf, output, options="", stranded="no",input_type="sam"):
    return Job(
        [input],
        [output],
        [['htseq_count', 'module_python']],
        command="""\
htseq-count {options} \\
  --stranded={stranded} \\
  --format={input_type} \\
  {input} \\
  {gtf}{output}""".format(
        options=options,
        stranded=stranded,
        input=input,
        gtf=gtf,
        output=" \\\n  > " + output if output else "",
        input_type=input_type
        )
    )
