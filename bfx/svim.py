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


def svim_ont(input_bam, output_directory):
    """
    Create a SVIM SV calling job for nanopore reads.

    :return: a job for nanopore SV calling with SVIM
    """

    genome_fasta = global_config_parser.param('DEFAULT', 'genome_fasta', required=True)

    min_mapq = global_config_parser.param('svim', 'min_map_qual')

    return Job(
        [input_bam],
        [os.path.join(output_directory, "final_results.vcf")],
        [["svim", "module_python3"]],
        command="""\
mkdir -p {output_directory} && \\
svim alignment --min_mapq {min_mapq} {other_options} {output_directory} {input_bam} {genome_fasta}
        """.format(
            output_directory=output_directory,
            min_mapq=min_mapq,
            other_options=global_config_parser.param('svim', 'other_options', required=False),
            input_bam=input_bam,
            genome_fasta=genome_fasta
        )
    )
