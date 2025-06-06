################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

from ..core.job import Job

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
