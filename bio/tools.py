#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

## function for awk tools ##



## function for python tools ## 
def py_addLengthRay (scaffolds_File, length_file, output):
    job = Job([scaffolds_File, length_file], [output], [['DEFAULT' , 'module_tools'], ['DEFAULT' , 'module_python']])
    
    job.command = \
"""python addLengthRay.py \\
  -s {scaFile} \\
  -l {lenFile}""".format(
        scaFile=scaffolds_File,
        lenFile=fastq
    )
    
    return job
def py_blastMatchSca (scaffolds_File, blast_file, output):
    job = Job([scaffolds_File, blast_file], [output], [['DEFAULT' , 'module_tools'], ['DEFAULT' , 'module_python']])
    
    job.command = \
"""python blastMatchSca.py \\
  -f {scaFile} \\
  -b {blastFile}""".format(
        scaFile=scaffolds_File,
        blastFile=blast_file
    )
    
    return job

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




## function for R tools ##

def r_select_scaffolds(input, output, name_sample, type_insert, min_insert_size=200):
    job = Job([input], [output], [['DEFAULT' , 'module_R']])
    
    job.command = \
""" R --no-save --args \\
  ./ \\
  ./ \\
  {name_sample} \\
  {type_insert} \\
  {min_insert_size} \\
  < analyse_selectSca.r \\""".format(
        name_sample=name_sample,
        type_insert=type_insert,
        min_insert_size=min_insert_size
    )
    
    return job

def r_find_cluster(input, output, unmap_type, name_sample, type_insert, max_insert_size=200, min_mapping_quality=10):
    job = Job([input], [output], [['DEFAULT' , 'module_R']])
    
    job.command = \
""" R --no-save --args \\
  ./ \\
  ./ \\
  {name_sample} \\
  {type_insert} \\
  {min_mapping_quality} \\
  {max_insert_size} \\
  < analyse_findCluster_{unmap_type}.r \\""".format(
        name_sample=name_sample,
        type_insert=type_insert,
        min_mapping_quality=min_mapping_quality,
        max_insert_size=max_insert_size,
        unmap_type=unmap_type
    )
    
    return job

def r_find_insert(input, output, name_sample, type_insert, mean_coverage=20, max_insert_size=200, min_overlap=2, exclu_file="None"):
    job = Job([input], [output], [['DEFAULT' , 'module_R']])
    
    job.command = \
""" R --no-save --args \\
  ./ \\
  ./ \\
  {name_sample} \\
  {type_insert} \\
  {mean_coverage} \\
  {max_insert_size} \\
  {min_overlap} \\
  {exclu_file} \\
  < analyse_findInsert.r \\""".format(
        name_sample=name_sample,
        type_insert=type_insert,
        mean_coverage=mean_coverage,
        max_insert_size=max_insert_size,
        min_overlap=min_overlap,
        exclu_file=exclu_file
    )
    
    return job

def r_filter_insert(input, output, name_sample, type_insert, mean_coverage=20, max_insert_size=200, strand=1, min_num_read=1, mean_read_length=100):
    job = Job([input], [output], [['DEFAULT' , 'module_R']])
    
    job.command = \
""" R --no-save --args \\
  ./ \\
  ./ \\
  {name_sample} \\
  {type_insert} \\
  {mean_coverage} \\
  {max_insert_size} \\
  {strand} \\
  {min_num_read} \\
  {mean_read_length} \\
  < analyse_findInsert.r \\""".format(
        name_sample=name_sample,
        type_insert=type_insert,
        mean_coverage=mean_coverage,
        max_insert_size=max_insert_size,
        strand=strand,
        min_num_read=min_num_read,
        mean_read_length=mean_read_length
    )
    
    return job



