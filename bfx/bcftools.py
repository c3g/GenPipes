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

def concat(inputs, output, options=None):
    """
    Concatenate or combine VCF/BCF files
    """
    if not isinstance(inputs, list):
        inputs=[inputs]
    
    return Job(
        inputs,
        [output],
        [
            ['bcftools_concat', 'module_bcftools']
        ],
        command="""\
bcftools \\
  concat -a {options} \\
  {inputs} \\
  {output}""".format(
        options=options if options else "",
        inputs="".join(" \\\n  " + input for input in inputs),
        output=" \\\n > " + output if output else ""
        )
    )

def view(input, output, filter_options):
    """
    Generalized view 
    """
    return Job(
        [input],
        [output],
        [
            ['bcftools_view', 'module_bcftools']
        ],
        command="""\
bcftools \\
  view {filter_options} \\
  {input}{output}""".format(
        filter_options=filter_options,
        input=" \\\n " + input if input else "",
        output=" \\\n > " + output if output else ""
        )
    )

def filter(input, output, filter_options):
    """
    Generalized filter function
    """
    return Job(
        [input],
        [output],
        [
            ['bcftools_filter', 'module_bcftools']
        ],
        command="""\
bcftools \\
  filter -m '+' -O v {filter_options} \\
  {input}{output}""".format(
        input=" \\\n " + input if input else "",
        filter_options=filter_options,
        output=" \\\n  > " + output if output else ""
        )
    )

