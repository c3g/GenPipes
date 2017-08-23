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


def bedGraphToBigWig(output_bed_graph, output_wiggle,header=True):
    
    if header :
        remove_head_command="""\
head -n 1  {output_bed_graph} > {output_bed_graph}.head.tmp && \\
awk ' NR > 1 ' {output_bed_graph} | sort -k1,1 -k2,2n > {output_bed_graph}.body.tmp && \\
cat {output_bed_graph}.head.tmp {output_bed_graph}.body.tmp > {output_bed_graph}.sorted && \\
rm {output_bed_graph}.head.tmp {output_bed_graph}.body.tmp""".format(
            output_bed_graph=output_bed_graph
        )
    else:
        remove_head_command="""\
sort -k1,1 -k2,2n {output_bed_graph} > {output_bed_graph}.sorted""".format(
            output_bed_graph=output_bed_graph
        )
    
    return Job(
        [output_bed_graph],
        [output_bed_graph+".sorted", output_wiggle],
        [
            ['bedtools', 'module_ucsc']
        ],
        command="""\
{remove_head_command} && \\
bedGraphToBigWig \\
  {output_bed_graph}.sorted \\
  {chromosome_size} \\
  {output_wiggle}""".format(
            remove_head_command=remove_head_command,
            chromosome_size=config.param('bedtools', 'chromosome_size', type='filepath'),
            output_bed_graph=output_bed_graph,
            output_wiggle=output_wiggle
        ),
        removable_files=[output_bed_graph+".sorted"]
    )

def bedToBigBed(bed_file, bigBed_file):
    
    return Job(
        [bed_file],
        [bigBed_file],
        [
            ['bedtools', 'module_ucsc']
        ],
        command="""\
bedToBigBed \\
  {bed_file} \\
  {chromosome_size} \\
  {bigBed_file}""".format(
            chromosome_size=config.param('bedtools', 'chromosome_size', type='filepath'),
            bed_file=bed_file,
            bigBed_file=bigBed_file
        ),
    )

