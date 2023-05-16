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

def trim( input1, input2, prefix, adapter_file, quality_offset):
    output_pair1 = prefix + "-trimmed-pair1.fastq.gz"
    output_pair2 = prefix + "-trimmed-pair2.fastq.gz"
    output_log = prefix + "-trimmed.log"

    if input2:  # Paired end reads
        inputs = [input1, input2]
        output = [output_pair1, output_pair2, output_log]
    else:   # Single end reads
        inputs = [input1]
        output = [output_pair1, output_log]

    return Job(
        inputs,
        output,
        [
            ['skewer_trimming', 'module_skewer']
        ],

        command="""\
$SKEWER_HOME/./skewer --threads {threads} {options} \\
  {adapter_file} \\
  {inputs} \\
  {outputs} \\
  -f {quality_offset}""".format(
        threads=config.param('skewer_trimming', 'threads', param_type='posint'),
        options=config.param('skewer_trimming', 'options'),
        adapter_file="-x " + adapter_file, 
        inputs=" \\\n  ".join(inputs),
        outputs="-o " + prefix,
        quality_offset="solexa" if quality_offset == 64 else "sanger",
        ),
    )
