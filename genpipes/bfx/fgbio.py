################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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

# Python Standard Modules
import logging

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

log = logging.getLogger(__name__)

def addumi(
    input_bam,
    input_umi,
    output_bam,
    output_bai,
    ini_section='fgbio_addumi'
    ):

    inputs = [input_bam, input_umi]
    outputs = [output_bam,output_bai]
    return Job(
        inputs,
        outputs,
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_fgbio']
        ],

        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $FGBIO_JAR AnnotateBamWithUmis \\
  --input {input_bam} \\
  --fastq {input_umi} \\
  --output {output_bam} \\
  {other_options}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            input_bam=input_bam,
            input_umi=input_umi,
            output_bam=output_bam,
            other_options=global_conf.global_get(ini_section, 'other_options', required=False)
            ),
        removable_files=[output_bam]
        )


def correct_readname(
    input_umi,
    output_umi_corrected
    ):
    
    inputs = [input_umi]
    outputs = [output_umi_corrected]
    
    if input_umi.lower().endswith('.gz'):
        input_opener="zcat "
    else:
        input_opener="cat "
        
    return Job(
        inputs,
        outputs,        
        command="""\
{input_opener} {input_umi} | tr ' ' '_' | gzip -c - > {output_umi_corrected}""".format(
        input_opener=input_opener,
        input_umi=input_umi,
        output_umi_corrected=output_umi_corrected
        ),
        removable_files=[output_umi_corrected]
    )


def trim_primers(input_bam, output_bam, hard_clip=False, ini_section='fgbio_trim_primers'):
    inputs = [input_bam]
    outputs = [output_bam]

    return Job(
        inputs,
        outputs,
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_fgbio']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $FGBIO_JAR TrimPrimers \\
  --input {input_bam} \\
  --primers {primers} \\
  --output {output_bam} \\
  {hard_clip} \\
  {other_options}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            input_bam=input_bam,
            primers=global_conf.global_get(ini_section, 'primers'),
            output_bam=output_bam,
            hard_clip="-H" if hard_clip else "",
            other_options=global_conf.global_get(ini_section, 'other_options')
            )
        )
