#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def compute_tdf(input, output):
    job = Job([input], [output], [['compute_tdf', 'moduleVersion.java'], ['compute_tdf', 'moduleVersion.igvtools']])

    job.command = \
"""igvtools count -f min,max,mean \\
  {input} \\
  {output} \\
  {genome}""".format(
        input=input,
        output=output,
        genome=config.param('compute_tdf', 'igvGenome')
    )

    return job
