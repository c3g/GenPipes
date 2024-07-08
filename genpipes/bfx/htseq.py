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
import os

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def htseq_count(input, gtf, output, options="", stranded="no",input_type="bam"):
    return Job(
        [input],
        [output],
        [['htseq_count', 'module_htseq']],
        command="""\
htseq-count {input} \\
  {options} \\
  --stranded={stranded} \\
  --format={input_type} \\
  {gtf}{output}""".format(
        options=options,
        stranded=stranded,
        input=input,
        gtf=gtf,
        output=" \\\n  > " + output if output else "",
        input_type="-" if input_type == "/dev/stdin" else input_type
        )
    )
