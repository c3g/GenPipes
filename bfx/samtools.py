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

# MUGQIC Modules
from core.config import *
from core.job import *

def index(input):
    return Job(
        [input],
        [input + ".bai"],
        [['samtools_index', 'module_samtools']],
        command="""\
samtools index \\
  {input}""".format(
        input=input
        )
    )

def faidx(input, filter=None):
    return Job(
        [input],
        [input + ".fai"],
        [['samtools_index', 'module_samtools']],
        command="""\
samtools faidx \\
  {input}{filter}""".format(
        input=input,
        filter=filter if filter else ""
        )
    )

def flagstat(input, output):
    return Job(
        [input],
        [output],
        [['samtools_flagstat', 'module_samtools']],
        command="""\
samtools flagstat \\
  {input} \\
  > {output}""".format(
        input=input,
        output=output
        ),
        removable_files=[output]
    )

def mpileup(input_bams, output, other_options="", region=None, regionFile=None, ini_section='rawmpileup'):

    return Job(
        input_bams,
        [output],
        [
            [ini_section, 'module_samtools']
        ],
        command="""\
samtools mpileup {other_options} \\
  -f {reference_fasta}{region}{regionFile}{input_bams}{output}""".format(
        other_options=other_options,
        reference_fasta=config.param('samtools_mpileup', 'genome_fasta', type='filepath'),
        regionFile=" \\\n  -l " + regionFile if regionFile else "",
        region=" \\\n  -r " + region if region else "",
        input_bams="".join([" \\\n  " + input_bam for input_bam in input_bams]),
        output=" \\\n  > " + output if output else ""
        )
    )



def merge(sample_output, input_bams):
    """
    merges an array of bams into a single bam
    """


    command = "samtools merge {sample_output} {input_bams}".format(sample_output = sample_output, input_bams = " ".join(map(str.strip, input_bams)))


    return Job( input_files = input_bams,
                    output_files = [sample_output],
                    module_entries = [['hicup_align', 'module_samtools']],
                    command = command
                )



def sort(input_bam, output_prefix, sort_by_name=False):
    output_bam = output_prefix + ".bam"

    return Job(
        [input_bam],
        [output_bam],
        [['samtools_sort', 'module_samtools']],
        command="""\
samtools sort {other_options}{sort_by_name} \\
  {input_bam} \\
  {output_prefix}""".format(
        other_options=config.param('samtools_sort', 'other_options', required=False),
        sort_by_name=" -n" if sort_by_name else "",
        input_bam=input_bam,
        output_prefix=output_prefix if config.param('samtools_sort', 'module_samtools').split("/")[2] == "0.1.19" else "> " + output_prefix + ".bam"
        ),
        removable_files=[output_bam]
    )

def view(input, output=None, options=""):
    return Job(
        [input],
        [output],
        [
            ['samtools_view', 'module_samtools']
        ],
        command="""\
samtools view {options} \\
  {input}{output}""".format(
        options=options,
        input=input,
        output=" \\\n  > " + output if output else ""
        ),
        removable_files=[output]
    )

def bcftools_cat(inputs, output):
    if not isinstance(inputs, list):
        inputs=[inputs]
    return Job(
        inputs,
        [output],
        [['bcftools_cat', 'module_bcftools']],
        command="""\
bcftools concat \\
  {inputs}{output}""".format(
        inputs=" \\\n  ".join(inputs),
        output=" \\\n  > " + output if output else ""
        )
    )

def bcftools_view(input, output, options="", pair_calling=False):
    return Job(
        [input],
        [output],
        [['bcftools_view', 'module_bcftools']],
        command="""\
bcftools view {pair_calling} {options} \\
  {input}{output}""".format(
        options=options,
        pair_calling="-T pair" if pair_calling else "",
        input=input,
        output=" \\\n  > " + output if output else ""
        )
    )

def bcftools_call(input, output, options="", pair_calling=False):
    return Job(
        [input],
        [output],
        [['bcftools_call', 'module_bcftools']],
        command="""\
bcftools call {pair_calling} {options} \\
  {input}{output}""".format(
        options=options,
        pair_calling="-T pair" if pair_calling else "",
        input=input,
        output=" \\\n  > " + output if output else ""
        )
    )
  
def bcftools_call_pair(input, output, options="", pair_calling=False):
    return Job(
        [input],
        [output],
        [['samtools_paired', 'module_samtools']],
        command="""\
$BCFTOOLS_BIN/bcftools view {pair_calling} {options} \\
  {input}{output}""".format(
        options=options,
        pair_calling="-T pair" if pair_calling else "",
        input=input,
        output=" \\\n  > " + output if output else ""
        )
    )

def bcftools_cat_pair(inputs, output):
    if not isinstance(inputs, list):
        inputs=[inputs]
    return Job(
        inputs,
        [output],
        [['samtools_paired', 'module_samtools']],
        command="""\
$BCFTOOLS_BIN/bcftools cat \\
  {inputs}{output}""".format(
        inputs=" \\\n  ".join(inputs),
        output=" \\\n  > " + output if output else ""
        )
    )

def bcftools_view_pair(input, output, options="", pair_calling=False):
    return Job(
        [input],
        [output],
        [['samtools_paired', 'module_samtools']],
        command="""\
$BCFTOOLS_BIN/bcftools view {pair_calling} {options} \\
  {input}{output}""".format(
        options=options,
        pair_calling="-T pair" if pair_calling else "",
        input=input,
        output=" \\\n  > " + output if output else ""
        )
    )
