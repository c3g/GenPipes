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


def salmon_index(reference_fasta,
                 output_directory):
    """
    Create a Salmon index from the reference fasta.
    """

    return Job(
        [reference_fasta],
        [output_directory],
        [["salmon_index", "module_salmon"]],
        command="""\
salmon index \\
 -t {reference_fasta} \\
 -i {index_directory}""".format(
            reference_fasta=reference_fasta,
            index_directory=output_directory),
        removable_files=[output_directory]
    )


def salmon_quant(readset_name,
                 read1_fastq, read2_fastq,
                 output_directory,
                 salmon_index):
    """
    Quantify transcripts using Salmon.
    """
    out_sf=os.path.join(output_directory,readset_name,readset_name, ".sf")

    return Job(
        [read1_fastq, read2_fastq],
        [out_sf],
        [["salmon_quant", "module_salmon"]],
        command="""\
salmon quant \\
 -i {salmon_index} \\
 -l A \\
 -o {output_directory} \\ 
 -1 {read1_fastq} \\
 -2 {read2_fastq} \\
 {other_options}""".format(
            salmon_index=salmon_index,
            output_directory=output_directory,
            read1_fastq=read1_fastq,
            read2_fastq=read2_fastq,
            other_options=config.param('salmon', 'other_options', required=False)
        ),
        removable_files=[]
    )
