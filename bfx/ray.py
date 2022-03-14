################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

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
        
    return Job(
        inputList, 
        [#TODO add output files
         os.path.join(pathOut, "Scaffolds.fasta"),
         os.path.join(pathOut, "ScaffoldLengths.txt")
        ], 
        [
            ['ray_assembly', 'module_gcc'],
            ['ray_assembly', 'module_openmpi'],
            ['ray_assembly', 'module_ray']
        ],
        command="""\
mpiexec {optionmpi} Ray -k {kmer} \\
  {pair} \\
  {single} \\
  -o {output}""".format(
        optionmpi=global_config_parser.param('ray_assembly', 'openmpi_options', required=False),
        rmdir=pathOut,
        mkdir=pathOut,
        kmer=global_config_parser.param('ray_assembly', 'kmer'),
        options=options,
        pair=strPair,
        single=strSingle,
        output=pathOut
        )
    )
