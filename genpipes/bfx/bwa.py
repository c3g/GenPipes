################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes. If not, see <http://www.gnu.org/licenses/>.
################################################################################

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

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

    other_options = global_conf.global_get(
        ini_section, 'bwa_other_options',
        required=False) if global_conf.global_get(
        ini_section, 'bwa_other_options', required=False) else ""

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
            idxbase=ref if ref else global_conf.global_get(ini_section, 'genome_bwa_index', param_type='filepath'),
            in1fastq=in1fastq,
            in2fastq=" \\\n  " + in2fastq if in2fastq else "",
            out_sam=" \\\n  > " + out_sam if out_sam else ""
        ),
        removable_files=[out_sam]
    )
