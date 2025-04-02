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


def align(
    input_file,
    read_group,
    sample_name,
    out_bam=None,
    sort=False,
    ini_section='pbmm2_align'
    ):
    """
    Align revio reads to a reference using pbmm2.

    :return: a job for revio alignment
    """

    genome_fasta = global_conf.global_get(ini_section, 'genome_fasta', required=True)

    return Job(
        [input_file],
        [out_bam, out_bam + ".bai"],
        [
            [ini_section, "module_pbmm2"]
        ],
        command="""\
pbmm2 align {other_options} \\
  --preset {pbmm2_preset} --log-level INFO \\
  -j {threads} -J {sort_threads} -m 4G {sort} \\
  --sample {sample_name} \\
  {read_group} \\
  {genome_fasta} \\
  {input_file} \\
  {out_bam}""".format(
            pbmm2_preset=global_conf.global_get(ini_section, 'preset'),
            threads=global_conf.global_get(ini_section, 'threads'),
            sort_threads=global_conf.global_get(ini_section, 'sort_threads'),
            sort="--sort " if sort else "",
            sample_name=sample_name,
            read_group="--rg " + read_group if read_group else "",
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            genome_fasta=genome_fasta,
            input_file=input_file,
            out_bam=out_bam if out_bam else ""
        ),
        removable_files=[out_bam]
    )
