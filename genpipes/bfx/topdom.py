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

def create_input(input_matrix, tmp_matrix, output_matrix, output_script, res):
    ## make TopDom R script:
    FileContent="""
source(\\\"{script}\\\"); \\
TopDom(matrix.file=\'{tmp_matrix}\', window.size={n}, outFile=\'{output_matrix}\')""".format(
        script=os.path.expandvars("${R_TOOLS}/TopDom_v0.0.2.R"),
        tmp_matrix=tmp_matrix,
        n=global_conf.global_get('identify_TADs', 'TopDom_n'),
        output_matrix=output_matrix
    )

    return Job(
        [input_matrix],
        [tmp_matrix],
        [
            ["identify_TADs", "module_R"],
            ["identify_TADs", "module_mugqic_tools"]
        ],
        command="""
echo \"{FileContent}\" > {fileName} && \\
{script} {input_matrix} {res}""".format(
            FileContent=FileContent,
            fileName=output_script,
            script="CreateTopDomMat.sh",
            input_matrix=input_matrix,
            res=res
        ),
        removable_files=[tmp_matrix]
    )

def call_TADs(input_matrix, matrix, R_script):
    return Job(
        [input_matrix],
        [matrix + ".bed", matrix + ".binSignal", matrix + ".domain"],
        [
            ["identify_TADs", "module_R"],
            ["identify_TADs", "module_mugqic_tools"]
        ],
        command="""
Rscript {fileName} && \\
rm {fileName}""".format(
            fileName=R_script
        ),
        removable_files=[input_matrix]
    )
