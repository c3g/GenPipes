################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import os

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def ray(pathOut, pair1=[], pair2=[], single=[], options=""):
    if len(pair1)!=len(pair2):
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
        optionmpi=global_conf.global_get('ray_assembly', 'openmpi_options', required=False),
        rmdir=pathOut,
        mkdir=pathOut,
        kmer=global_conf.global_get('ray_assembly', 'kmer'),
        options=options,
        pair=strPair,
        single=strSingle,
        output=pathOut
        )
    )
