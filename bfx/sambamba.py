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

# MUGQIC Modules
from core.config import *
from core.job import *

def sort(input_bam,
         output_bam,
         tmp_dir,
         other_options=None):
    
    return Job(
        [input_bam],
        [output_bam],
        [
            ['sambamba_sort_sam', 'module_sambamba']
        ],
        command="""\
sambamba sort {options} \\
  {input} \\
  --tmpdir {tmp} \\
  {output}""".format(
        options=other_options if other_options else "",
        input=input_bam,
        output="--out " + output_bam if output_bam else "",
        tmp=tmp_dir
        )
    )

def index(input,
          output,
          other_options=None
          ):

    return Job(
        [input],
        [output],
        [
            ['sambamba_index', 'module_sambamba']
        ],
        command="""\
sambamba index {options} \\
  {input} \\
  {output}""".format(
        options=other_options if other_options else "",
        input=input,
        output=output,
        )
    )

def merge(input_bams,
          output_bam,
          ini_section='sambamba_merge_sam_files'
          ):

    return Job(
        input_bams,
        [output_bam],
        [
            [ini_section, 'module_samtools'],
            [ini_section, 'module_sambamba']
        ],
        command="""\
sambamba merge {options} \\
  {output} \\
  {input}""".format(
        options=global_config_parser.param(ini_section, 'options'),
        input="".join([" \\\n  " + input_bam for input_bam in input_bams]),
        output=output_bam
        )
    )

def markdup(input_bam,
            output_bam,
            tmp_dir,
            other_options=None
            ):
    if not isinstance(input_bam, list):
        input_bam=[input_bam]

    return Job(
        input_bam,
        [output_bam],
        [
            ['sambamba_mark_duplicates', 'module_samtools'],
            ['sambamba_mark_duplicates', 'module_sambamba']
        ],
        command="""\
sambamba markdup {other_options} \\
  {input} \\
  --tmpdir {tmp} \\
  {output}""".format(
        other_options=other_options,
        tmp=tmp_dir,
        input=" \\\n  ".join(input for input in input_bam),
        output=output_bam,
        )
    )

def view(input_bam,
         output_bam,
         options,
         chr=[]
         ):

    return Job(
        [input_bam],
        [output_bam],
        [
            ['DEFAULT', 'module_sambamba'],
        ],
        command="""\
sambamba view {options} \\
  {input} \\
  {output} {chr}""".format(
        options=options,
        input=input_bam,
        output="-o " + output_bam if output_bam else "",
        chr=chr if chr else "",
        )
    )

def flagstat(input, output, options=None):

    return Job(
        [input],
        [output],
        [
            ['DEFAULT', 'module_sambamba'],
        ],
        command="""\
sambamba flagstat {options} \\
  {input} \\
  {output}""".format(
        options=options if options else "",
        input=input,
        output="> " + output,
        )
    )
