################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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


def report(
    input,
    output_dir,
    tumor_id
    ):

    assembly = config.param('report_cpsr', 'assembly')
    output = [
        os.path.join(output_dir, tumor_id + ".cpsr." + assembly + ".json.gz"),
        os.path.join(output_dir, tumor_id + ".cpsr." + assembly + ".html")
    ]

    if config.param('report_cpsr', 'module_pcgr').split("/")[2] >= "1":
        call = 'cpsr'
        module = 'module_pcgr'
    else:
        call = 'cpsr.py'
        module = 'module_cpsr'

    return Job(
        [input],
        output,
        [
            ['report_cpsr', module],
        ],
        command="""\
{call} {options} \\
    --input_vcf {input} \\
    --pcgr_dir $PCGR_DATA \\
    --output_dir {output_dir} \\
    --genome_assembly {assembly} \\
    --sample_id {tumor_id}""".format(
            call=call,
            options=config.param('report_cpsr', 'options'),
            input=input,
            output_dir=output_dir,
            assembly=config.param('report_cpsr', 'assembly'),
            tumor_id=tumor_id
        )
    )

def parse_cpsr_passed_variants_pt(input_file):
    """
    Parse PCGR passed variants.
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export cpsr_passed_variants=`grep "cpsr-gene-annotate - INFO - Number of PASSed variant calls:" {input_file} | awk -F': ' '{{print $2}}'`"""
        )
