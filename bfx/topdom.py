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

def create_input(input_matrix, tmp_matrix, output_matrix, output_script, res):
    ## make TopDom R script:
    FileContent="""
source(\\\"{script}\\\"); \\
TopDom(matrix.file=\'{tmp_matrix}\', window.size={n}, outFile=\'{output_matrix}\')""".format(
        script=os.path.expandvars("${R_TOOLS}/TopDom_v0.0.2.R"),
        tmp_matrix=tmp_matrix,
        n=global_config_parser.param('identify_TADs', 'TopDom_n'),
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
