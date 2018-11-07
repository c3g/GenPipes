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
  --output {output_bam}
  {other_options}""".format(
        tmp_dir=config.param('fgbio_addumi', 'tmp_dir'),
        java_other_options=config.param('fgbio_addumi', 'java_other_options'),
        ram=config.param('fgbio_addumi', 'ram'),
        input_bam=input_bam
        input_umi=input_umi
        output_bam=output_bam
        ),
        removable_files=[output_bam]
    )
