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
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def ballgown( 
    input_files,
    design_file,
    output_dir
    ):

    return  Job(
        input_files + [design_file],
        [os.path.join(output_dir, "gene_exp.diff"), os.path.join(output_dir, "transcript_exp.diff")],
        [
            ['ballgown', 'module_mugqic_tools'],
            ['ballgown', 'module_R']
        ],
        command="""\
Rscript $R_TOOLS/ballgown.R \\
  -d {design_file} \\
  -o {output_dir} \\
  """.format(
        design_file=design_file,
        output_dir=output_dir,
    ))

