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
import logging
import os

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

def addumi(
    input_bam,
    input_umi,
    output_bam,
    output_bai
    ):

    inputs = [input_bam, input_umi]
    outputs = [output_bam,output_bai]
    return Job(
        inputs,
        outputs,
        [
            ['fgbio_addumi', 'module_java'],
            ['fgbio_addumi', 'module_fgbio']
        ],

        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $FGBIO_JAR AnnotateBamWithUmis \\
  --input {input_bam} \\
  --fastq {input_umi} \\
  --output {output_bam} \\
  {other_options}""".format(
            tmp_dir=global_config_parser.param('fgbio_addumi', 'tmp_dir'),
            java_other_options=global_config_parser.param('fgbio_addumi', 'java_other_options'),
            ram=global_config_parser.param('fgbio_addumi', 'ram'),
            input_bam=input_bam,
            input_umi=input_umi,
            output_bam=output_bam,
            other_options=global_config_parser.param('fgbio_addumi', 'other_options') if global_config_parser.param('fgbio_addumi', 'other_options', required=False) else ""
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


def trim_primers(input_bam, output_bam, hard_clip=False):
    inputs = [input_bam]
    outputs = [output_bam]

    return Job(
        inputs,
        outputs,
        [
            ['fgbio_trim_primers', 'module_java'],
            ['fgbio_trim_primers', 'module_fgbio']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $FGBIO_JAR TrimPrimers \\
  --input {input_bam} \\
  --primers {primers} \\
  --output {output_bam} \\
  {hard_clip} \\
  {other_options}""".format(
            tmp_dir=global_config_parser.param('fgbio_trim_primers', 'tmp_dir'),
            java_other_options=global_config_parser.param('fgbio_trim_primers', 'java_other_options'),
            ram=global_config_parser.param('fgbio_trim_primers', 'ram'),
            input_bam=input_bam,
            primers=global_config_parser.param('fgbio_trim_primers', 'primers'),
            output_bam=output_bam,
            hard_clip="-H" if hard_clip else "",
            other_options=global_config_parser.param('fgbio_trim_primers', 'other_options')
            )
        )
