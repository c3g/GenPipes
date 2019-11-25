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

def calculate_reproducible_score(
    output_dir,
    sample1,
    sample2,
    chr
    ):

    return Job(
        [sample1,sample2],
        [os.path.join("deseq_results.csv")],
        [
            ['differential_expression_deseq', 'module_mugqic_tools'],
            ['differential_expression_deseq', 'module_R']
        ],
        command="""\
Rscript /home/pubudu/projects/def-bourqueg/pubudu/job_outputs/hicrep.R \\
  -1 {sample1} \\
  -2 {sample2} \\
  -c {chr} \\
  -o {output_dir}""".format(
        sample1=sample1,
        sample2=sample2,
        output_dir=output_dir,
        chr=chr
    ))
