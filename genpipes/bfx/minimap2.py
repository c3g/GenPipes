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


def minimap2_ont(
    input,
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
        [input],
        [out_sam],
        [
            [ini_section, "module_minimap2"]
        ],
        command="""\
minimap2 \\
  -t {threads} \\
  -ax {minimap_preset} {other_options} \\
  -R {read_group} \\
  {genome_fasta} \\
  {input}{out_sam}""".format(
            threads=global_conf.global_get(ini_section, 'threads'),
            minimap_preset=global_conf.global_get(ini_section, 'preset'),
            read_group=read_group,
            other_options=global_conf.global_get(ini_section, 'minimap2_other_options', required=False),
            genome_fasta=genome_fasta,
            input=input + "/*.fastq*" if input else "-",
            out_sam=" \\\n  > " + out_sam if out_sam else ""
        ),
        removable_files=[out_sam]
    )
