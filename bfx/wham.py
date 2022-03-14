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

def call_sv(input_normal, input_tumor, output):
    if input_tumor is not None:
        inputs = [input_normal, input_tumor]
        
    else:
        inputs = [input_normal]

    return Job(
        inputs,
        [output],
        [
            ['wham_call_sv', 'module_wham']
        ],
        command="""\
export EXCLUDE={exclude}
WHAM-GRAPHENING \\
    -k -z   \\
    -x {cores}  \\
    -e $EXCLUDE \\
    -a {genome} \\
    -f {input_tumor}{input_normal} \\
    {output}""".format(
            exclude=global_config_parser.param('wham_call_sv', 'exclude'),
            cores=global_config_parser.param('wham_call_sv', 'cores'),
            genome=global_config_parser.param('wham_call_sv', 'genome_fasta', param_type='filepath'),
            input_tumor=input_tumor + "," if input_tumor else "",
            input_normal=input_normal,
            output=" \\\n  > " + output if output else ""
        )
    )

def merge(input_vcf, output_vcf):
    return Job(
        [input_vcf],
        [output_vcf],
        [
            ['wham_call_sv', 'module_wham']
        ],
        command="""\
mergeIndvs  \\
    -f {input_vcf}      \\
    {output_vcf}""".format(
            input_vcf=input_vcf,
            output_vcf=" \\\n > " + output_vcf if output_vcf else ""
        )
    )


def genotype(input_vcf, input_normal, input_tumor, output):
    if input_tumor is not None:
        inputs = [input_vcf, input_normal, input_tumor]
    
    else:
        inputs = [input_vcf, input_normal]
        
    return Job(
        inputs,
        [output],
        [
            ['wham_call_sv', 'module_wham']
        ],
        command="""\
WHAM-GRAPHENING \\
    -b {input_vcf}  \\
    -x {cores}  \\
    -a {genome} \\
    -f {input_tumor}{input_normal} \\
    {output}""".format(
            input_vcf=input_vcf,
            cores=global_config_parser.param('wham_call_sv', 'cores'),
            genome=global_config_parser.param('wham_call_sv', 'genome_fasta', param_type='filepath'),
            input_tumor=input_tumor + "," if input_tumor else "",
            input_normal=input_normal,
            output=" \\\n  > " + output if output else ""
        )
    )
