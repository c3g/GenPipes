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

def sym_link(input, readset, out_dir, type=None):
    sample = ""
    if type == "raw_reads":
        sample = readset.sample.name

    else:
        sample = readset.name


    prefix = os.path.join(out_dir, "deliverables", sample, global_config_parser.param('DEFAULT', 'experiment_type_abrev'), type)
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
        prefix=prefix,
        input=os.path.join(out_dir,input),
        output=output
        )
    )

def sym_link_pair(input, tumor_pair, out_dir, type=None, sample=None, profyle=False):
    if profyle:
        pair = ""
        if not (type == "raw_reads" or type == "alignment"):
            pair = tumor_pair.pair_profyle + "/"

        if sample == "Normal":
            prefix = os.path.join(out_dir, "analyses", tumor_pair.name, tumor_pair.normal_profyle, global_config_parser.param('DEFAULT', 'experiment_type_abrev'), pair + type)

        else:
            prefix = os.path.join(out_dir, "analyses", tumor_pair.name, tumor_pair.tumor_profyle, global_config_parser.param('DEFAULT', 'experiment_type_abrev'), pair + type)

    else:
        if sample == "Normal":
            prefix = os.path.join(out_dir, "deliverables", tumor_pair.name, tumor_pair.normal.name, global_config_parser.param('DEFAULT', 'experiment_type_abrev'), type)

        else:
            prefix = os.path.join(out_dir, "deliverables", tumor_pair.name, tumor_pair.tumor.name, global_config_parser.param('DEFAULT', 'experiment_type_abrev'), type)

    input_postfix = input.split("/")[-1]
    output = os.path.join(prefix, input_postfix)

    return Job(
        [input],
        [output],
        command="""\
mkdir -p {prefix} && \\
ln -s -f \\
  {input} \\
  {output}""".format(
        prefix=prefix,
        input=os.path.join(out_dir, input),
        output=output
        )
    )

def md5sum(input, output, out_dir):
    return Job(
        [input],
        [output],
    command="""\
md5sum {input} \\
  > {output}""".format(
        input=os.path.join(out_dir, input),
        output=os.path.join(out_dir, output),
        )
    )
