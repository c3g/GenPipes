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

# MUGQIC Modules
from core.config import *
from core.job import *

def add_reject(input, output):
    """
    Adds REJECT to the filter. Used with scalpel when merging common indels
    """
    return Job(
        [input],
        [output],
        [
            ['bcftools_add_reject', 'module_bcftools']
        ],
        command="""\
bcftools \\
  filter -m '+' -O v --soft-filter 'REJECT' -e '%TYPE="indel"' \\
  {input}{output}""".format(
        input=input,
        output=" \\\n  > " + output if output else ""
        )
    )

def add_chi2Filter(input, output):
    """
    Adds CHI2FILTER to the filter. Used with scalpel on somatic indels
    """
    return Job(
        [input],
        [output],
        [
            ['bcftools_add_chi2Filter', 'module_bcftools']
        ],
        command="""\
bcftools \\
  filter -m '+' -O v --soft-filter 'CHI2FILTER' -e 'INFO/CHI2 > 20.0' \\
  {input}{output}""".format(
        input=input,
        output=" \\\n  > " + output if output else ""
        )
    )


