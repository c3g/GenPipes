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

def decompose_and_normalize_mnps(inputs, vt_output=None):
    if not isinstance(inputs, list):
        inputs = [inputs]

    return Job(
        inputs,
        [vt_output],
        [
            ['decompose_and_normalize_mnps', 'module_htslib'],
            ['decompose_and_normalize_mnps', 'module_vt']
        ],
        command="""\
zless {input} | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r {reference_sequence} - | sed -e '#0/.#0/0#' \\
        {vt_output}""".format(
        input=" \\\n  ".join(input for input in inputs),
        reference_sequence=global_config_parser.param('decompose_and_normalize_mnps', 'genome_fasta', param_type='filepath'),
        vt_output="> " + vt_output if vt_output else " ",
        )
    )

def sort(input, output, options):

    return Job(
        [input],
        [output],
        [
            ['vt_sort', 'module_htslib'],
            ['vt_sort', 'module_vt']
        ],
        command="""\
vt sort {options} -o {output} {input}""".format(
        options=options,
        input=" \\\n " + input if input else "-",
        output=output
        )
    )
