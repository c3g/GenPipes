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


def report(input, output_dir, tumor_id):
    assembly = config.param('cpsr_report', 'assembly')
    output_json = os.path.join(output_dir, tumor_id + ".cpsr." + assembly + ".json.gz"),
    output_html = os.path.join(output_dir, tumor_id + ".cpsr." + assembly + "html"),
    
    return Job(
        input,
        [output_json, output_html],
        [
            ['cpsr_report', 'mugqic_cpsr'],
        ],
        command="""\
cpsr.py {options} \\
    --input_vcf {input} \\
    --pcgr_dir $PCGR_DATA \\
    --output_dir {output_dir} \\
    --genome_assembly {assembly} \\
    --sample_id {tumor_id}""".format(
            options=config.param('cpsr_report', 'option'),
            input=input,
            output_dir=output_dir,
            assembly=config.param('cpsr_report', 'assembly'),
            tumor_id=tumor_id
        )
    )