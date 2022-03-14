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

from core.config import *
from core.job import *

def call(input_pair, output, sv_type, genotype_file=None):
    return Job(
        input_pair,
        [output],
        [
            ['delly_call_filter', 'module_delly']
        ],
        command="""\
export OMP_NUM_THREADS={threads} && export LC_ALL=C && \\
delly call {options} \\
    -t {sv_type}    \\
    {exclude_list}  \\
    -o {output}     \\
    -g {genome}     \\
    {genotype_file} \\
    {input_pair}""".format(
            options=global_config_parser.param('delly_call_filter', 'options'),
            threads=global_config_parser.param('delly_call_filter', 'threads'),
            sv_type=sv_type,
            exclude_list="\\\n    -x " + global_config_parser.param('delly_call_filter', 'exclude_list') if global_config_parser.param('delly_call_filter', 'exclude_list') else "",
            output=output,
            genome=global_config_parser.param('delly_call_filter', 'genome_fasta', param_type='filepath'),
            genotype_file="-v " + genotype_file if genotype_file else "",
            input_pair="".join(" \\\n " + sample for sample in input_pair)
        )
    )

def filter(input_bcf, output_bcf, sv_type, type, options, sample_file):
    if not isinstance(input_bcf, list):
        input_bcf=[input_bcf]
    return Job(
        input_bcf,
        [output_bcf],
        [
            ['delly_call_filter', 'module_delly']
        ],
        command="""\
delly filter {options} \\
    -t {sv_type}            \\
    -f {type} {output_bcf}  \\
    {sample_file}           \\
    {input_bcf}""".format(
            options=options if options else "",
            sv_type=sv_type,
            type=type,
            output_bcf="-o " + output_bcf if output_bcf else "",
            sample_file="-s " + sample_file if sample_file else "",
            input_bcf=input_bcf
        )
    )

def merge(input_bcfs, output, sv_type):
    if not isinstance(input_bcfs, list):
        input_bcfs=[input_bcfs]
    return Job(
        input_bcfs,
        [output],
        [
            ['delly_call_filter', 'module_delly']
        ],
        command="""\
delly merge {options} \\
    -t {sv_type}    \\
    -o {output}     \\
    {input_bcfs}""".format(
            options=global_config_parser.param('delly_merge_germline', 'options') if global_config_parser.param('delly_merge_germline', 'options') else "",
            sv_type=sv_type,
            output=output,
            input_bcfs="".join(" \\\n " + bcf for bcf in input_bcfs)
        )
    )
