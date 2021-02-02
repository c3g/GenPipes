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
import logging

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)


def freebayes(input_bam, output_file, options=None, ini_section='freebayes'):
    # inputs = [input]
    # output = [prefix + ".vcf"]

    return Job(
        input_files=[input_bam],
        output_files=[output_file],
        module_entries=[
            [ini_section, 'module_freebayes']
        ],

        command="""\
freebayes -f {reference_genome} \\
  {bed_targets} \\
  {options} \\
  {other_options} \\
  {input_bam} \\
  > {output_file}""".format(
    reference_genome=config.param(ini_section, 'genome_fasta', type='filepath'),
    bed_targets="-t " + config.param(ini_section, 'bed_targets', type='filepath') if config.param(ini_section, 'bed_targets') else "",
    options=options,
    other_options=config.param(ini_section, 'other_options'),
    input_bam=input_bam,
    output_file=output_file
    ),
        )
