#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
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

def hicrep(
    base_dir,
    sample1,
    sample2,
    chr
    ):

    return  Job(
        [count_matrix],
        [os.path.join(output_dir, "deseq_results.csv"), os.path.join(output_dir, "dge_results.csv")],
        [
            ['differential_expression_deseq', 'module_mugqic_tools'],
            ['differential_expression_deseq', 'module_R']
        ],
        command="""\
Rscript $R_TOOLS/deseq.R \\
  -d {design_file} \\
  -c {count_matrix} \\
  -o {output_dir} \\
  {localfit}""".format(
        design_file=design_file,
        count_matrix=count_matrix,
        output_dir=output_dir,
        localfit=localfit
    ))

def edger(
    design_file,
    count_matrix,
    output_dir
    ):

    return  Job(
        [count_matrix],
        [os.path.join(output_dir, "edger_results.csv")],
        [
            ['differential_expression_edger', 'module_mugqic_tools'],
            ['differential_expression_edger', 'module_R']
        ],
        command="""\
Rscript $R_TOOLS/edger.R \\
  -d {design_file} \\
  -c {count_matrix} \\
  -o {output_dir}""".format(
        design_file=design_file,
        count_matrix=count_matrix,
        output_dir=output_dir
    ))