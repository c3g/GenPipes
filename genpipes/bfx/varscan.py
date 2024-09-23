################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
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
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################


# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def mpileupcns(
        input,
        output,
        sample_list,
        ini_section='germline_varscan2'
):

    return Job(
        [input, sample_list],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_varscan'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $VARSCAN2_JAR mpileup2cns {other_options} \\
  {input} \\
  --output-vcf 1 \\
  --vcf-sample-list {sampleNames}{output}""".format(
        tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
        java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
        ram=global_conf.global_get(ini_section, 'ram'),
        other_options=global_conf.global_get(ini_section, 'mpileup_other_options'),
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
            [ini_section, 'module_java'],
            [ini_section, 'module_varscan'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $VARSCAN2_JAR somatic \\
  {input_pair} \\
  {output} \\
  {other_options} \\
  --output-vcf 1 --mpileup 1""".format(
        tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
        java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
        ram=global_conf.global_get(ini_section, 'ram'),
        other_options=global_conf.global_get(ini_section, 'other_options'),
        input_pair=input_pair,
        output=output,
        )
    )

def fpfilter_somatic(
        input_vcf,
        input_readcount,
        output=None,
        ini_section='varscan2_somatic'
):

    return Job(
        [input_vcf, input_readcount],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_varscan'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $VARSCAN2_JAR fpfilter \\
  {input_vcf} \\
  {input_readcount} \\
  {options} \\
  {output}""".format(
        tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
        java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
        ram=global_conf.global_get(ini_section, 'ram'),
        options=global_conf.global_get(ini_section, 'fpfilter_options'),
        input_vcf=input_vcf,
        input_readcount=input_readcount,
        output="--output-file " + output if output else ""
        )
    )
