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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job


def run(
    input_bam,
    output_stats,
    output_blocks,
    output_summary,
    input_vcf1,
    output_vcf1,
    input_vcf2=None,
    output_vcf2=None,
    input_vcf3=None,
    output_vcf3=None,
    ini_section='hiphase'
    ):
    """
    Joint phasing of small, structural, and tandem repeat variants
    for PacBio sequencing data with HiPhase.

    :return: a job for HiPhase phasing
    """

    genome_fasta = global_conf.global_get(ini_section, 'genome_fasta', required=True)
    outputs = [
        output_stats,
        output_blocks,
        output_summary,
        output_vcf1,
        output_vcf2,
        output_vcf3
    ]

    return Job(
        [input_bam, input_vcf1, input_vcf2, input_vcf3],
        outputs,
        [
            [ini_section, "module_hiphase"]
        ],
        command="""\
hiphase {other_options} \\
  --reference {genome_fasta} \\
  --bam {input_file} \\
  --threads {threads} \\
  --stats-file {output_stats} \\
  --blocks-file {output_blocks} \\
  --summary-file {output_summary} \\
  --vcf {input_vcf1} \\
  --output-vcf {output_vcf1} \\
  {input_vcf2} \\
  {output_vcf2} \\
  {input_vcf3} \\
  {output_vcf3}""".format(
            threads=global_conf.global_get(ini_section, 'threads'),
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            genome_fasta=genome_fasta,
            input_file=input_bam,
            output_stats=output_stats,
            output_blocks=output_blocks,
            output_summary=output_summary,
            input_vcf1=input_vcf1,
            output_vcf1=output_vcf1,
            input_vcf2="--vcf " + input_vcf2 if input_vcf2 else "",
            output_vcf2="--output-vcf " + output_vcf2 if output_vcf2 else "",
            input_vcf3="--vcf " + input_vcf3 if input_vcf3 else "",
            output_vcf3="--output-vcf " + output_vcf3 if output_vcf3 else ""
        ),
    )
