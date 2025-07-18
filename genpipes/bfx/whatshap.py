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


def haplotag(
    input_bam,
    input_vcf,
    output_bam,
    ini_section='whatshap'
    ):
    """
    Create a haplotagged bam file with Whatshap.

    :return: a job for whatshap haplotagging
    """

    return Job(
        [input_bam, input_vcf],
        [output_bam],
        [
            [ini_section, "module_whatshap"],
        ],
        command="""\
whatshap haplotag \\
  -o {output_bam} \\
  {other_options} \\
  --reference {genome_fasta} \\
  {input_vcf} \\
  {input_bam}""".format(
            input_bam=input_bam,
            input_vcf=input_vcf,
            output_bam=output_bam,
            genome_fasta=global_conf.global_get(ini_section, 'genome_fasta'),
            other_options=global_conf.global_get(ini_section, 'other_options', required=False)
        )
    )
