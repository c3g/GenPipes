#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def ray(pathOut, pair1=[], pair2=[], single=[], options=""):
    if (len(pair1)!=len(pair2)):
      raise Exception("Error Ray: the number of paired files is different :" + len(pair1) + " - " + len(pair2))
    
    inputList = []
    strPair=""
    for i in range(len(pair1)):
        strPair = strPair + "-p " + pair1[i] + " " + pair2[i] + " " 
        inputList.append(pair1[i])
        inputList.append(pair2[i])
    
    strSingle=""
    for i in range(len(single)):
        strSingle = strSingle + "-s " + single[i] + " "
        inputList.append(single[i])
        
    job = Job(
        inputList, 
        [#TODO add output files
         os.path.join(pathOut, "Scaffolds.fasta"),
         os.path.join(pathOut, "ScaffoldLengths.txt")
        ], 
        [
         ['DEFAULT', 'module_gcc'],
         ['DEFAULT', 'module_openmpi'],
         ['DEFAULT', 'module_ray']
        ]
    )
    
    job.command = \
"""rm -r {rmdir}; mpiexec {optionmpi} Ray -k {kmer} \\
  {pair} \\
  {single} \\
  -o {output}""".format(
        optionmpi=config.param('DEFAULT', 'openmpi_options', required=False),
        rmdir=pathOut,
        mkdir=pathOut,
        kmer=config.param('assembly_of_unmap', 'kmer'),
        options=options,
        pair=strPair,
        single=strSingle,
        output=pathOut
    )

    return job
