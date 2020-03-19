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
import re

# MUGQIC Modules
from core.config import config
from core.job import Job, concat_jobs

from bfx import bash_cmd as bash


def fastqc(
    input1,
    input2,
    output,
    adapter_file=None,
    use_tmp=False
    ):

    if input2:  # Paired end reads
        inputs = [input1, input2]
    else:       # Single end reads
        inputs = [input1]

    outputs = [output]

    output_directory = os.path.dirname(output)
    tmp_directory = output_directory + ".tmp"

    (input_basename, file_format) = os.path.splitext(input1)
    file_format = re.sub("^\.", "", file_format)
    if file_format == 'gz':
        (input_basename, file_format) = os.path.splitext(input_basename)
        file_format = re.sub("^\.", "", file_format)

    return Job(
        inputs,
        outputs,
        [
            ['fastqc', 'module_fastqc'],
            ['fastqc', 'module_java']
        ],
        command="""\
fastqc \\
  --outdir {output_directory} \\
  --threads {threads} \\
  {adapter} \\
  --format {file_format} \\
  {tmp} \\
  {inputs} && \\
rm -r {tmp}""".format(
            output_directory=output_directory,
            threads=config.param('fastqc', 'threads', param_type='posint'),
            adapter="--adapters " + adapter_file if adapter_file else "",
            file_format=file_format,
            tmp="-d " + tmp_directory if tmp_directory else "",
            inputs="\\\n  ".join(inputs)
        ),
        removable_files=[]
    )
