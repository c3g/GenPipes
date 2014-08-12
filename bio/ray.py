#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def ray(pathOut, pair1=[], pair2=[], single=[], options=""):
    if (len(pair1)!=len(pair2)):
      raise Exception("Error Ray: the number of paired files is different :" + len(pair1) + " - " + len(pair2))
    
    inputList = []
    pair=""
    for i in range(len(pair1)):
        pair = pair + "-p " + pair1[i] + " " + pair2[i] + " " 
        inputList.append(pair1[i])
        inputList.append(pair2[i])
    
    single=""
    for element in single:
        single = single + "-s " + element + " "
        inputList.append(element)
        
    job = Job(
        inputList, 
        [#TODO add output files
         os.path.join(pathOut, "Scaffolds.fasta"),
         os.path.join(pathOut, "Scaffolds.fasta.length")
        ], 
        [
         ['DEFAULT', 'module_gcc'],
         ['DEFAULT', 'module_openmpi'],
         ['DEFAULT', 'module_ray']
        ]
    )
    
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
