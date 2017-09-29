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

def bedGraphToBigWig(input_bed_graph, output_wiggle, header=True):

    # Check it the input is a real bedGrah (i.e. contains the bedGraph header : track type=bedGraph)
    # or if it is just a regular bed file (i.e. no bedGraph header)
    if header :
        if os.path.splitext(input_bed_graph)[1] == ".gz":
            remove_head_command="""\
zcat {input_bed_graph} | head -n 1 > {input_bed_graph}.head.tmp && \\
zcat {input_bed_graph} | awk ' NR > 1 ' | sort  --temporary-directory={temp_dir} -k1,1 -k2,2n > {input_bed_graph}.body.tmp && \\
cat {input_bed_graph}.head.tmp {input_bed_graph}.body.tmp > {input_bed_graph}.sorted && \\
rm {input_bed_graph}.head.tmp {input_bed_graph}.body.tmp""".format(
                input_bed_graph=input_bed_graph,
                temp_dir=config.param('ucsc', 'tmp_dir',  required=True)
            )
        else:
            remove_head_command="""\
head -n 1  {input_bed_graph} > {input_bed_graph}.head.tmp && \\
awk ' NR > 1 ' {input_bed_graph} | sort  --temporary-directory={temp_dir} -k1,1 -k2,2n > {input_bed_graph}.body.tmp && \\
cat {input_bed_graph}.head.tmp {input_bed_graph}.body.tmp > {input_bed_graph}.sorted && \\
rm {input_bed_graph}.head.tmp {input_bed_graph}.body.tmp""".format(
                input_bed_graph=input_bed_graph,
                temp_dir=config.param('ucsc', 'tmp_dir',  required=True)
            )
    else:
        if os.path.splitext(input_bed_graph)[1] == ".gz":
          remove_head_command="""\
zcat {input_bed_graph} | sort --temporary-directory={temp_dir} -k1,1 -k2,2n > {input_bed_graph}.sorted""".format(
                input_bed_graph=input_bed_graph,
                temp_dir=config.param('ucsc', 'tmp_dir',  required=True)
            )
        else:
            remove_head_command="""\
sort --temporary-directory={temp_dir} -k1,1 -k2,2n {input_bed_graph} > {input_bed_graph}.sorted""".format(
                input_bed_graph=input_bed_graph,
                temp_dir=config.param('ucsc', 'tmp_dir', required=True)
            )

    return Job(
        [input_bed_graph],
        [input_bed_graph+".sorted", output_wiggle],
        [
            ['ucsc', 'module_ucsc']
        ],
        command="""\
{remove_head_command} && \\
bedGraphToBigWig \\
  {input_bed_graph}.sorted \\
  {chromosome_size} \\
  {output_wiggle}""".format(
            remove_head_command=remove_head_command,
            chromosome_size=config.param('ucsc', 'chromosome_size', type='filepath'),
            input_bed_graph=input_bed_graph,
            output_wiggle=output_wiggle
        ),
        removable_files=[input_bed_graph+".sorted"]
    )

def bedToBigBed(bed_file, bigBed_file):

    return Job(
        [bed_file],
        [bigBed_file],
        [
            ['ucsc', 'module_ucsc']
        ],
        command="""\
bedToBigBed \\
  {bed_file} \\
  {chromosome_size} \\
  {bigBed_file}""".format(
            chromosome_size=config.param('ucsc', 'chromosome_size', type='filepath'),
            bed_file=bed_file,
            bigBed_file=bigBed_file
        )
    )
