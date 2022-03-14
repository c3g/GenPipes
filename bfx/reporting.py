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
from core.config import global_config_parser
from core.job import Job
from utils import utils

def formatvcf(input_vcf, output_vcf):
    return Job(
        [input_vcf],
        [output_vcf],
        [
	        ['format_vcf', 'module_mugqic_tools'],
            ['format_vcf', 'module_python']
        ],
        command="""\
python $PYTHONTOOLS/format2pcgr.py {options} \\
        -i {input}  \\
        -o {output}""".format(
            input=input_vcf,
			output=output_vcf,
            options=global_config_parser.param('format_vcf', 'options'),
        )
    )
