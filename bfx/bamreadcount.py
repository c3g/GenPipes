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

def readcount(input, bed, output):

    return Job(
        [input],
        [output],
        [
            [ 'varscan2_readcount_fpfilter', 'module_bamreadcount'],
        ],
        command="""\
$BAMREADCOUNT_BIN/bam-readcount {options} \\
  -f {reference_sequence} \\
  {input} \\
  -l {bed} \\
  {output}""".format(
        options=global_config_parser.param('varscan2_readcount_fpfilter', 'readcount_options'),
        reference_sequence=global_config_parser.param('varscan2_readcount_fpfilter', 'genome_fasta', param_type='filepath'),
        input=input,
        bed=bed,
        output="> " + output if output else ""
        )
    )

