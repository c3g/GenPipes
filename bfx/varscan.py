################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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

def mpileupcns(
        input,
        output,
        sample_list
):

    return Job(
        [input, sample_list],
        [output],
        [
            ['germline_varscan2', 'module_java'],
            ['germline_varscan2', 'module_varscan'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $VARSCAN2_JAR mpileup2cns {other_options} \\
  {input} \\
  --output-vcf 1 \\
  --vcf-sample-list {sampleNames}{output}""".format(
        tmp_dir=config.param('germline_varscan2', 'tmp_dir'),
        java_other_options=config.param('germline_varscan2', 'java_other_options'),
        ram=config.param('germline_varscan2', 'ram'),
        other_options=config.param('germline_varscan2', 'mpileup_other_options'),
        input=" \\\n " + input if input else "",
        sampleNames=sample_list,
        output=" \\\n  > " + output if output else ""
        )
    )

def somatic(
        input_pair,
        output,
        output_vcf_dep=[],
        output_snp_dep=[],
        output_indel_dep=[],
        ini_section='varscan2_somatic'
):

    return Job(
        [input_pair],
        [output_vcf_dep, output_snp_dep, output_indel_dep],
        [
            ['varscan2_somatic', 'module_java'],
            ['varscan2_somatic', 'module_varscan'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $VARSCAN2_JAR somatic \\
  {input_pair} \\
  {output} \\
  {other_options} \\
  --output-vcf 1 --mpileup 1""".format(
        tmp_dir=config.param(ini_section, 'tmp_dir'),
        java_other_options=config.param(ini_section, 'java_other_options'),
        ram=config.param(ini_section, 'ram'),
        other_options=config.param(ini_section, 'other_options'),
        input_pair=input_pair,
        output=output,
        )
    )

def fpfilter_somatic(input_vcf, input_readcount, output=None):

    return Job(
        [input_vcf, input_readcount],
        [output],
        [
            ['varscan2_somatic', 'module_java'],
            ['varscan2_somatic', 'module_varscan'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $VARSCAN2_JAR fpfilter \\
  {input_vcf} \\
  {input_readcount} \\
  {options} \\
  {output}""".format(
        tmp_dir=config.param('varscan2_somatic', 'tmp_dir'),
        java_other_options=config.param('varscan2_somatic', 'java_other_options'),
        ram=config.param('varscan2_somatic', 'ram'),
        options=config.param('varscan2_somatic', 'fpfilter_options'),
        input_vcf=input_vcf,
        input_readcount=input_readcount,
        output="--output-file " + output if output else ""
        )
    )
