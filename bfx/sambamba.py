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

#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def sort(input_bam, output_bam, tmp_dir):

    return Job(
        [input_bam],
        [output_bam],
        [
            ['sambamba_sort_sam', 'module_samtools'],
            ['sambamba_sort_sam', 'module_sambamba']
        ],
        command="""\
sambamba sort {options} \\
  {input} \\
  --tmpdir {temp} \\
  {output}""".format(
        options=config.param('sambamba_sort_sam', 'options'),
        input=input_bam,
        output="--out " + output_bam if output_bam else "",
        temp=tmp_dir
        )
    )

def index(input_bam, tmp_dir, output_dep=[]):

    return Job(
        [input_bam],
        [output_dep],
        [
            ['sambamba_index', 'module_samtools'],
            ['sambamba_index', 'module_sambamba']
        ],
        command="""\
sambamba index {options} \\
  {input}""".format(
        options=config.param('sambamba_index', 'options'),
        input=input,
        )
    )

def merge(input_bams, output_bam):

    return Job(
        input_bams,
        [output_bam],
        [
            ['sambamba_merge_sam_files', 'module_samtools'],
            ['sambamba_merge_sam_files', 'module_sambamba']
        ],
        command="""\
sambamba merge {options} \\
  {output} \\
  {input}""".format(
        options=config.param('sambamba_merge_sam_files', 'options'),
        input="".join([" \\\n  " + input_bam for input_bam in input_bams]),
        output=output_bam
        )
    )

def markdup(input_bam, output_bam, tmp_dir):

    return Job(
        [input_bam],
        [output_bam],
        [
            ['sambamba_mark_duplicates', 'module_samtools'],
            ['sambamba_mark_duplicates', 'module_sambamba']
        ],
        command="""\
sambamba markdup {options} \\
  {input} \\
  --tmpdir {temp} \\
  {output}""".format(
        options=config.param('sambamba_mark_duplicates', 'options'),
        input=input_bam,
        output=output_bam,
        temp=tmp_dir
        )
    )
