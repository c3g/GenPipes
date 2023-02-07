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
import os
import re

from core.job import *

def bigWigToBedGraph(input_bigWigFile, output_bedgraph):
    # Remove the chr_ prefix to convert on whole genome or change chrom number to convert only for specified chromosome
    # bigWigFile_name = os.path.basename(bigWigFile)
    # output_bedgraph = os.path.join(output_dir, "chr1_"+bigWigFile_name+".bedgraph")

    return Job(
        [input_bigWigFile],
        [],
        [['ucsc', 'module_ucsc']],
        name="bigwig_to_bedgraph",
        command="""bigWigToBedGraph {bigwig} {output_file}""".format(
            bigwig=input_bigWigFile,
            output_file=output_bedgraph
            )
        )

def bigWigInfo(input_bigwig, output):
   # output = os.path.join(output_dir, "bigwiginfo_"+ os.path.basename(input_bigwig) + ".txt")

    return Job(
        [input_bigwig],
        [output],
        [['ucsc', 'module_ucsc']],
        name="bigwiginfo",
        command="""bigWigInfo {input_bigwig} > {output}""".format(
            input_bigwig=input_bigwig,
            output=output
            )
        )