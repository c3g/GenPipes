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


def minimap2_ont(readset_name,
                 read_fastq_dir,
                 output_directory):
    """
    Align nanopore reads to a reference using minimap2.

    :return: a job for nanopore alignment
    """

    genome_fasta = config.param('minimap2_align', 'genome_fasta', required=True)

    minimap_preset = config.param('minimap2_align', 'preset')

    out_sam = os.path.join(output_directory, "Aligned.out.sam")
    out_bam = os.path.join(output_directory, readset_name + ".sorted.bam")
    out_bai = os.path.join(output_directory, readset_name + ".sorted.bai")

    return Job(
        [read_fastq_dir],
        [out_bam, out_bai, out_sam],
        [["minimap2", "module_minimap2"],
         ["minimap2", "module_samtools"]],
        command="""\
mkdir -p {output_directory} && \\
minimap2 -ax {minimap_preset} {other_options} {genome_fasta} {read_fastq_dir}/*.fastq > {out_sam} && \\
samtools view -b {out_sam} -@ 8 | samtools sort -@ 6 --reference {genome_fasta} -T {output_directory} > {out_bam} && \\ 
samtools index {out_bam} {out_bai}
        """.format(
            output_directory=output_directory,
            minimap_preset=minimap_preset,
            other_options=config.param('minimap2_align', 'other_options', required=False),
            genome_fasta=genome_fasta,
            read_fastq_dir=read_fastq_dir,
            out_sam=out_sam,
            out_bam=out_bam,
            out_bai=out_bai
        ),
        removable_files=[out_sam]
    )
