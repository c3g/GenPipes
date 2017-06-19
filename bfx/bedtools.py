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
import os

# MUGQIC Modules
from core.config import *
from core.job import *

from bfx import ucsc

def graph(input_bam, output_bed_graph, output_wiggle, library_type="PAIRED_END"):

    if library_type == "PAIRED_END":
        if "forward" in output_bed_graph:
            samtools_options = "-F 256 -f 81 "
        elif "reverse" in output_bed_graph:
            samtools_options = "-F 256 -f 97 "
        else:
            raise Exception("Error: PAIRED_END library was provided but no strand orientation could be determined from " + output_bed_graph + "...")
    else:
        samtools_options = "-F 256"

    return concat_jobs([
        Job(
            [input_bam],
            [output_bed_graph, output_wiggle],
            [
                ['bedtools_graph', 'module_samtools'],
                ['bedtools_graph', 'module_bedtools']
            ],
            command="""\
nmblines=$(samtools view {samtools_options} {input_bam} | wc -l) && \\
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \\
genomeCoverageBed -bg -split -scale $scalefactor \\
  -ibam {input_bam} \\
  -g {chromosome_size} \\
  > {output_bed_graph}""".format(
                samtools_options=samtools_options,
                input_bam=input_bam,
                chromosome_size=config.param('bedtools_graph', 'chromosome_size', type='filepath'),
                other_options=config.param('bedtools_graph', 'other_options'),
                output_bed_graph=output_bed_graph
            )
        ),
        ucsc.bedgraph_to_bigbwig(
            output_bed_graph,
            output_wiggle
        )
    ])

def intersect(input_bam, output_bam, target_bed):

    return Job(
        [input_bam],
        [output_bam],
        [
            ['bedtools_intersect', 'module_bedtools']
        ],
        command="""\
bedtools intersect \\
  -a {input_bam} \\
  -b {target_bed} \\
  {other_options} > {output_bam}""".format(
            input_bam=input_bam,
            target_bed=target_bed,
            other_options=config.param('bedtools_intersect', 'other_options'),
            output_bam=output_bam
        )
    )

def bamtobed(input_bam, output_bed):

    return Job(
        [input_bam],
        [output_bed],
        [
            ['bedtools_intersect', 'module_bedtools']
        ],
        command="""\
bedtools bamtobed \\
  -i {input_bam} > {output_bed}""".format(
            input_bam=input_bam,
            output_bed=output_bed
        )
    )

def coverage(output_bed, output_file):

    return Job(
        [output_bed],
        [output_file],
        [
            ['bedtools_intersect', 'module_bedtools']
        ],
        command="""\
bedtools coverage \\
  -b {output_bed} \\
  -a {intervals} \\
  {other_options} > {output_file}""".format(
            output_bed=output_bed,
            intervals=config.param('bedtools_coverage', 'gc_intervals'),
            other_options=config.param('bedtools_coverage', 'other_options'),
            output_file=output_file
        )
    )