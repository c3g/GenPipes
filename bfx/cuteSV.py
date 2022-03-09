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

# MUGQIC Modules
from core.config import *
from core.job import *


def cuteSV_ont(input_bam, sample_name, output_directory):
    """
    Create a cuteSV SV calling job for nanopore reads.

    :return: a job for nanopore SV calling with cuteSV
    """

    genome_fasta = config.param('DEFAULT', 'genome_fasta', required=True)
    out_vcf = os.path.join(output_directory, "cuteSV_structural_variants.vcf")

    return Job(
        [input_bam],
        [os.path.join(output_directory, "cuteSV_structural_variants.vcf")],
        [["cuteSV", "module_cuteSV"]],
        command="""\
mkdir -p {output_directory} && \\
cuteSV {other_options} -t {threads} \\
    --max_cluster_bias_INS {max_cluster_bias_ins} \\
    --diff_ratio_merging_INS {diff_ratio_merging_ins} \\
    --max_cluster_bias_DEL {max_cluster_bias_del} \\
    --diff_ratio_merging_DEL {diff_ratio_merging_del} \\
    --report_readid \\
    --genotype \\
    -S {sample_name} \\
    {input_bam} \\
    {genome_fasta} \\
    {out_vcf} \\
    {tmp_dir} 
        """.format(
            output_directory=output_directory,
            other_options=config.param('cuteSV', 'other_options', required=False),
            threads=config.param('cuteSV', 'threads', required=True),
            max_cluster_bias_ins=config.param('cuteSV', 'max_cluster_bias_ins', required=True),
            diff_ratio_merging_ins=config.param('cuteSV', 'diff_ratio_merging_ins', required=True),
            max_cluster_bias_del=config.param('cuteSV', 'max_cluster_bias_del', required=True),
            diff_ratio_merging_del=config.param('cuteSV', 'diff_ratio_merging_del', required=True),
            sample_name=sample_name,
            input_bam=input_bam,
            genome_fasta=genome_fasta,
            out_vcf=out_vcf,
            tmp_dir=config.param('DEFAULT', 'tmp_dir', required=True)
        )
    )
