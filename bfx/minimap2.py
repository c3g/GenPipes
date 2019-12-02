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


def minimap2_ont(read_fastq_dir,
                 output_directory,
                 sort_bam=True):
    """
    Align nanopore reads to a reference using minimap2.

    :return: a job for nanopore alignment
    """

    genome_fasta = config.param('minimap2', 'genome_fasta', required=True)

    minimap_preset = config.param('minimap2', 'preset')

    bam_name = "Aligned.sortedByCoord.out.bam"
    out_bam = os.path.join(output_directory, bam_name)

    return Job(
        [read_fastq_dir],
        [out_bam],
        [["minimap2", "module_minimap2"]],
        command="""\
mkdir -p {output_directory} && \\
cd {output_directory} && \\
minimap2 -ax {minimap_preset} {other_options} {genome_fasta} {read_fastq_dir}/*.fastq > Aligned.out.sam && \\
samtools view -b Aligned.out.sam -@ 8 | samtools sort -@ 6 --reference {genome_fasta} > {bam_name} && \\ 
samtools index {bam_name}
        """.format(
            output_directory=output_directory,
            minimap_preset=minimap_preset,
            other_options=config.param('minimap2', 'other_options', required=False),
            genome_fasta=genome_fasta,
            read_fastq_dir=read_fastq_dir,
            bam_name=bam_name
        )
    )
