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

import os

from core.job import *

def bigWigToBedGraph(input_dir, bigWigFile, output_dir):
    # Remove the chr_ prefix to convert on whole genome or change chrom number to convert only for specified chromosome
    bigWigFile_name = os.path.basename(bigWigFile)
    output_bedgraph = os.path.join(output_dir, "chr1_"+bigWigFile_name+".bedgraph")

    return Job(
        [input_dir],
        [output_dir],
        [['ucsc', 'module_ucsc']],
        name = "bigwig_to_bedgraph",
        command = """\
bigWigToBedGraph \\
  {bigwig} \\
  {output_file} 2> {stderr}; [ -s {stderr} ] || rm -f {stderr}""".format( # If there are no errors, we delete the stderr file
        bigwig = bigWigFile,
        output_file = output_bedgraph,
        stderr = os.path.join(output_dir+"_error", bigWigFile_name+".error")
        )
    )

def bigWigInfo(input_dir, bigWigFile, output_dir):
    output = os.path.join(output_dir, "bigwiginfo_"+ os.path.basename(bigWigFile) + ".txt")

    return Job(
        [input_dir],
        [output],
        [['ucsc', 'module_ucsc']],
        name = "bigwiginfo",
        command = """\
bigWigInfo \\
  {bigWigFile} &> {output}""".format(
        bigWigFile = bigWigFile,
        output = output
        )
    )