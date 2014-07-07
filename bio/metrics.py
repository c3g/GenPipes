#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def dna_sample_metrics(input_directory, output, experiment_type="unknown"):
    job = Job([input_directory], [output], [['dna_sample_metrics', 'moduleVersion.R'], ['dna_sample_metrics', 'moduleVersion.tools']])

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

def snv_graph_metrics(list, output_basename):
    job = Job([list], [output_basename + ".snvGraphMetrics_listFiles.txt"], [['snv_graph_metrics', 'moduleVersion.R'], ['snv_graph_metrics', 'moduleVersion.tools']])

    job.command = \
"""Rscript \$R_TOOLS/snvGraphMetrics.R \\
  {list} \\
  {output_basename}""".format(
        list=list,
        output_basename=output_basename
    )

    return job

def vcf_stats(input, output, list):
    job = Job([input, list], [output], [['vcf_stats', 'moduleVersion.python'], ['vcf_stats', 'moduleVersion.tools']])

    job.command = \
"""python \$PYTHON_TOOLS/vcfStats.py \\
  -v {input} \\
  -d {dictionary} \\
  -o {output} \\
  -f {list}""".format(
        input=input,
        dictionary=config.param('vcf_stats', 'referenceSequenceDictionary', type='filepath'),
        output=output,
        list=list
    )

    return job
