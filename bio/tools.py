#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

## function for awk tools ##
def awk_samToFastq (input=None, output):
    cmd = "sh samToFastq.awk " + input + " > " + output if input else "sh samToFastq.awk - > " + output
    job = Job(output_files=[output], command=cmd, module_entries=[['DEFAULT' , 'module_tools']])

    return job
    

## function for python tools ## 
def py_equalFastqFile (fastq_ref, fastq, output):
    job = Job([fastq_ref, fastq], [output], [['DEFAULT' , 'module_tools'], ['DEFAULT' , 'module_python']])
    
    job.command = \
"""python equalFastqFile.py \\
  -r {ref} \\
  -f {fastq}""".format(
        ref=fastq_ref,
        fastq=fastq
    )

    return job

## function for R tools ##


## function for perl tools ##

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
