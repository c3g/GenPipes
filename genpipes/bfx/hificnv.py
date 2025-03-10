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
    input_file,
    output_prefix,
    sample_name,
    input_maf = None,
    ini_section='hificnv'
    ):
    """
    Copy number variant calling and visualisation with HiFiCNV.

    :return: a job for sawfish discover
    """

    genome_fasta = global_conf.global_get(ini_section, 'genome_fasta', required=True)
    outputs = [
        f"{output_prefix}.{sample_name}.copynum.bedgraph",
        f"{output_prefix}.{sample_name}.depth.bw",
        f"{output_prefix}.{sample_name}.maf.bw",
        f"{output_prefix}.{sample_name}.vcf.gz"
    ]

    return Job(
        [input_file, input_maf],
        outputs,
        [
            [ini_section, "module_hificnv"]
        ],
        command="""\
hificnv {other_options} \\
  --bam {input_file} \\
  --ref {genome_fasta} \\
  {exclude} \\
  {expected_cn} \\
  {input_maf} \\
  --threads {threads} \\
  --output-prefix {output_prefix}""".format(
            threads=global_conf.global_get(ini_section, 'threads'),
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            genome_fasta=genome_fasta,
            input_maf="--maf " + input_maf if input_maf else "",
            exclude="--exclude " + global_conf.global_get(ini_section, 'excluded_regions', required=False) if global_conf.global_get(ini_section, 'excluded_regions', required=False) else "",
            expected_cn="--expected-cn " + global_conf.global_get(ini_section, 'cn_bed', required=False) if global_conf.global_get(ini_section, 'cn_bed', required=False) else "",
            input_file=input_file,
            output_prefix=output_prefix
        ),
    )
