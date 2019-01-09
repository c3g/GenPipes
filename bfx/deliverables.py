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

def sym_link(input, readset, type=None):
    sample = ""
    if type == "raw_reads":
        sample = readset.sample.name

    else:
        sample = readset.name


    prefix = os.path.join("deliverables", sample, config.param('DEFAULT','experiment_type_abrev'), type)
    input_postfix = input.split("/")[-1]

    output = os.path.join(prefix, input_postfix)

    return Job(
    [input],
    [output],
    command="""\
mkdir -p {prefix} && \\       
ln -sf \\
  {input} \\
  {output}""".format(
        prefix=os.path.abspath(prefix),
        input=os.path.abspath(input),
        output=os.path.abspath(output)
        )
    )

def sym_link_pair(input, tumor_pair, type=None, sample=None, profyle=False):
    
    if profyle == True:
        pair = ""
        if not (type == "raw_reads" or type == "alignment"):
            pair = tumor_pair.pair_profyle + "/"

        if sample == "Normal":
            prefix = os.path.join("analyses", tumor_pair.name, tumor_pair.normal_profyle, config.param('DEFAULT','experiment_type_abrev'), pair + type)

        else:
            prefix = os.path.join("analyses", tumor_pair.name, tumor_pair.tumor_profyle, config.param('DEFAULT','experiment_type_abrev'), pair + type)

    else:
        if sample == "Normal":
            prefix = os.path.join("deliverables", tumor_pair.name, tumor_pair.normal.name, config.param('DEFAULT','experiment_type_abrev'), type)

        else:
            prefix = os.path.join("deliverables", tumor_pair.name, tumor_pair.tumor.name, config.param('DEFAULT','experiment_type_abrev'), type)

    
    input_postfix = input.split("/")[-1]
    output = os.path.join(prefix, input_postfix)

    return Job(
        [input],
        [output],
        command="""\
mkdir -p {prefix} && \\       
ln -sf \\
  {input} \\
  {output}""".format(
        prefix=os.path.abspath(prefix),
        input=os.path.abspath(input),
        output=os.path.abspath(output)
        )
    )


def md5sum(input, output):
    return Job(
    [input],
    [output],
    command="""\
md5sum {input} > {output}""".format(
        input=os.path.abspath(input),
        output=os.path.abspath(output)
        )
    )
