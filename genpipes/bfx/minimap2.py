################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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

# Python Standard Modules
import os

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job


def minimap2_ont(
    read_fastq_dir,
    read_group,
    out_sam=None,
    ini_section='minimap2_ont'
    ):
    """
    Align nanopore reads to a reference using minimap2.

    :return: a job for nanopore alignment
    """

    genome_fasta = global_conf.global_get(ini_section, 'genome_fasta', required=True)

    return Job(
        [read_fastq_dir],
        [out_sam],
        [[ini_section, "module_minimap2"]],
        command="""\
minimap2 \\
  -ax {minimap_preset} \\
  -R {read_group} \\
  {other_options} \\
  {genome_fasta} \\
  {read_fastq_dir}/*.fastq*{out_sam}""".format(
            minimap_preset=global_conf.global_get(ini_section, 'preset'),
            read_group=read_group,
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            genome_fasta=genome_fasta,
            read_fastq_dir=read_fastq_dir,
            out_sam=" \\\n  > " + out_sam if out_sam else ""
        ),
        removable_files=[out_sam]
    )
