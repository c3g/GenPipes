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

def index(input):
    return Job(
        [input],
        [input + ".bwt"],
        [['bwa_mem', 'module_bwa']],
        command="""\
bwa index \\
  {input}""".format(
      input=input
      )
        )

def mem(in1fastq, in2fastq=None, out_sam=None, read_group=None, ref=None, ini_section='bwa_mem'):

    other_options = global_config_parser.param(
        ini_section, 'bwa_other_options',
        required=False) if global_config_parser.param(
        ini_section, 'bwa_other_options', required=False) else "",

    return Job(
        [in1fastq, in2fastq, ref + ".bwt" if ref else None],
        [out_sam],
        [[ini_section, "module_bwa"]],
        command="""\
bwa mem {other_options} \\
  {read_group} \\
  {idxbase} \\
  {in1fastq}{in2fastq}{out_sam}""".format(
      other_options=other_options,
      read_group=" \\\n  -R " + read_group if read_group else "",
      idxbase=ref if ref else global_config_parser.param(ini_section, 'genome_bwa_index', param_type='filepath'),
      in1fastq=in1fastq,
      in2fastq=" \\\n  " + in2fastq if in2fastq else "",
      out_sam=" \\\n  > " + out_sam if out_sam else ""
      ),
        removable_files=[out_sam]
        )
