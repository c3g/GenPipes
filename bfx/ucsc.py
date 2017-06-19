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

def bedgraph_to_bigbwig(input_bed_graph, output_wiggle, graph_header=False):

    if os.path.splitext(input_bed_graph)[1] == ".gz":
        sorted_bed_graph = re.sub(".gz", ".sorted", input_bed_graph)
        sort_command = """zcat {input_bed_graph} | sort -k1,1 -k2,2n > {sorted_bed_graph}""".format(
            input_bed_graph=input_bed_graph,
            sorted_bed_graph=sorted_bed_graph
        )
        if graph_header:
            sort_command = """zcat {input_bed_graph} | tail -n +2 | sort -k1,1 -k2,2n > {sorted_bed_graph}""".format(
                input_bed_graph=input_bed_graph,
                sorted_bed_graph=sorted_bed_graph
            )
    else:
        sorted_bed_graph = input_bed_graph + ".sorted"
        sort_command = """sort -k1,1 -k2,2n {input_bed_graph} > {sorted_bed_graph}""".format(
            input_bed_graph=input_bed_graph,
            sorted_bed_graph=sorted_bed_graph
        )
        if graph_header:
            sort_command = """cat {input_bed_graph} | tail -n +2 | sort -k1,1 -k2,2n > {sorted_bed_graph}""".format(
                input_bed_graph=input_bed_graph,
                sorted_bed_graph=sorted_bed_graph
            )

    return Job(
        [input_bed_graph],
        [output_wiggle],
        [
            ['ucsc', 'module_ucsc']
        ],
        command="""\
{sort_command} && \\
bedGraphToBigWig \\
  {sorted_bed_graph} \\
  {chromosome_size} \\
  {output_wiggle}""".format(
            sort_command=sort_command,
            chromosome_size=config.param('ucsc', 'chromosome_size', type='filepath'),
            sorted_bed_graph=sorted_bed_graph,
            output_wiggle=output_wiggle
        )
    )
