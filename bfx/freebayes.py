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


def freebayes_varcall(input_bams, output=None, region=None):

    return Job(
        input_bams,
        [output],
        [
            ['freebayes', 'module_freebayes'],
            ['freebayes', 'module_samtools']
        ],
        command="""\
freebayes -f {reference_fasta} -r {region} {freebayes_options} \\
{input_bams} > {output}  \\
""".format(
        tmp_dir=config.param('freebayes', 'tmp_dir'),
        reference_fasta=config.param('freebayes','genome_fasta',type='filepath'),
        region=" \\\n  " + region if region else "",
        freebayes_options=config.param("freebayes","freebayes_other_options"),
        input_bams="".join([" \\\n  " + input_bam for input_bam in input_bams]),
        output=" \\\n  " + output if output else ""
        )
    )


