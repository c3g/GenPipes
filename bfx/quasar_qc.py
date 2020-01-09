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

def create_fend_object(chromosome_lengths, output_file, res_chr):

    return Job(
        [chromosome_lengths],
        [output_file],
        [
            ['quality_scores', 'module_python']
        ],
        command="""\
            hifive fends \\
            -L {chromosome_lengths} \\
            --binned={res_chr} \\
            {output_file}""".format(
            chromosome_lengths=chromosome_lengths,
            output_file=output_file,
            res_chr=res_chr
        ))