#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def ray(pathOut, pair1=[], pair2=[], single=[], options=""):
    if (len(pair1)!=len(pair2)):
      raise Exception("Error Ray: the number of paired files is different :" + len(pair1) + " - " + len(pair2))
    
    input=[]
    pair=""
    for i in range(0, len(pair1)):
        pair = pair + "-p " + pair1[x] + " " + pair2[x]
        input.append(pair1[x])
        input.append(pair2[x])
    
    single=""
    for element in single:
	single = single + "-s " + element
        input.append(element)
        
    job = Job([input], [pathOut+"Scaffolds.fasta"], [['DEFAULT', 'module_gcc'],['DEFAULT', 'module_openmpi'],['DEFAULT', 'module_ray']])
    
    job.command = \
"""mpiexec Ray -k {kmer} \\
  {pair} \\
  {single} \\
  -o {output}""".format(
        kmer=config.param('ray', 'kmer'),
        options=options,
        pair=pair,
        single=single,
        output=pathOut + " \\\n  > "
    )

    return job
