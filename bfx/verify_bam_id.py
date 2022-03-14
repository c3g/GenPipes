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

def verify(input_bam, output_prefix):
    return Job(
        [input_bam],
        [output_prefix + ".selfSM"],
        [
            ['verify_bam_id', 'module_verify_bam_id']
        ],
        command="""\
verifyBamID \\
  --vcf {input_vcf} \\
  --bam {input_bam} \\
  --out {output_prefix} \\
  {other_options}""".format(
            input_vcf=global_config_parser.param('verify_bam_id', 'vcf', param_type='filepath'),
            input_bam=input_bam,
            output_prefix=output_prefix,
            other_options=global_config_parser.param('verify_bam_id', 'options')
        )
    )
