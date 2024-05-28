################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
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

# GenPipes Modules
from ..core.config import global_conf
from ..core.job import Job

def svim_ont(input_bam, output_directory):
    """
    Create a SVIM SV calling job for nanopore reads.

    :return: a job for nanopore SV calling with SVIM
    """

    genome_fasta = global_conf.global_get('DEFAULT', 'genome_fasta', required=True)

    min_mapq = global_conf.global_get('svim', 'min_map_qual')

    return Job(
        [input_bam],
        [os.path.join(output_directory, "final_results.vcf")],
        [["svim", "module_python"]],
        command="""\
mkdir -p {output_directory} && \\
svim alignment --min_mapq {min_mapq} {other_options} {output_directory} {input_bam} {genome_fasta}
        """.format(
            output_directory=output_directory,
            min_mapq=min_mapq,
            other_options=global_conf.global_get('svim', 'other_options', required=False),
            input_bam=input_bam,
            genome_fasta=genome_fasta
        )
    )
