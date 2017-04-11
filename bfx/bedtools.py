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

def graph(input_bam, output_bed_graph, output_wiggle, library_type="PAIRED_END"):

    if library_type == "PAIRED_END":
        samtools_options="-F 256 -f 81 "
    else :
        samtools_options="-F 256"

    return Job(
        [input_bam],
        [output_bed_graph, output_wiggle],
        [
            ['bedtools', 'module_samtools'],
            ['bedtools', 'module_bedtools'],
            ['bedtools', 'module_ucsc']
        ],
        command="""\
nmblines=$(samtools view {samtools_options} {input_bam} | wc -l) && \\
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \\
genomeCoverageBed -bg -split -scale $scalefactor \\
  -ibam {input_bam} \\
  -g {chromosome_size} \\
  > {output_bed_graph} && \\
sort -k1,1 -k2,2n {output_bed_graph} > {output_bed_graph}.sorted && \\
bedGraphToBigWig \\
  {output_bed_graph}.sorted \\
  {chromosome_size} \\
  {output_wiggle}""".format(
            samtools_options=samtools_options,
            input_bam=input_bam,
            chromosome_size=config.param('bedtools', 'chromosome_size', type='filepath'),
            output_bed_graph=output_bed_graph,
            output_wiggle=output_wiggle
        )
    )

def intersect(input_bam, output_bam, readset):
    target_bed = config.param('bedtools_intersect', 'target_bed', required=False)

    if target_bed:
        if target_bed == "auto":
            if readset.beds:
                target_bed = os.path.abspath(readset.beds[0])
        else:
            # Add filepath validation
            target_bed = config.param('bedtools_intersect', 'target_bed', type='filepath')

    return Job(
        [input_bam],
        [output_bam],
        [
            ['bedtools', 'module_bedtools']
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

def genomecov():

    if library_type == "PAIRED_END":
        samtools_options="-F 256 -f 81 "
    else :
        samtools_options="-F 256"

    return Job(
        [input_bam],
        [output_bed_graph, output_wiggle],
        [
            ['bedtools', 'module_samtools'],
            ['bedtools', 'module_bedtools']
        ],
        command="""\

  > {output_bed_graph} && \\
sort -k1,1 -k2,2n {output_bed_graph} > {output_bed_graph}.sorted && \\
bedGraphToBigWig \\
  {output_bed_graph}.sorted \\
  {chromosome_size} \\
  {output_wiggle}""".format(
            samtools_options=samtools_options,
            input_bam=input_bam,
            chromosome_size=config.param('bedtools', 'chromosome_size', type='filepath'),
            output_bed_graph=output_bed_graph,
            output_wiggle=output_wiggle
        )
    )
