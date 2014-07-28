#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def bed2interval_list(dictionary, bed, output):
    job = Job([dictionary, bed], [output], [['DEFAULT' , 'module_tools'], ['DEFAULT' , 'module_perl']])

    if not dictionary:
        dictionary = config.param('DEFAULT', 'genome_dictionary', type='filepath')

    job.command = \
"""bed2IntervalList.pl \\
  --dict {dictionary} \\
  --bed {bed} \\
  > {output}""".format(
        dictionary=dictionary,
        bed=bed,
        output=output
    )

    return job

def filter_long_indel(input, output):
    job = Job([input], [output], [['DEFAULT' , 'module_tools'], ['DEFAULT' , 'module_perl']])

    job.command = \
"""filterLongIndel.pl \\
  {input} \\
  > {output}""".format(
        input=input,
        output=output
    )

    return job
