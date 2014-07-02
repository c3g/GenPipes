#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def dna_sample_metrics(input_directory, output, experiment_type="unknown"):
    job = Job([input_directory], [output], [['metrics', 'moduleVersion.R'], ['metrics', 'moduleVersion.tools']])

    job.command = \
"""Rscript \$R_TOOLS/DNAsampleMetrics.R \\
  {input_directory} \\
  {output} \\
  {experiment_type}""".format(
        input_directory=input_directory,
        output=output,
        experiment_type=experiment_type
    )

    return job
