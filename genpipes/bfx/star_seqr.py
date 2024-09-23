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
import logging
import os

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

log = logging.getLogger(__name__)

def run(fastqs1, fastqs2, output_dir, ini_section='run_star_seqr'):
    if not isinstance(fastqs1, list):
        fastqs1 = [fastqs1]
        
    if not isinstance(fastqs2, list):
        fastqs2 = [fastqs2]
        
    output_file = os.path.join(output_dir + "_STAR-SEQR", "star_seqr_STAR-SEQR_candidates.txt")
    return Job(
        fastqs1,
        [output_file],
        [
            [ini_section, 'module_conda']
        ],

        command="""\
    starseqr.py -t {threads} {options} \\
      -i {genome_build} \\
      -g {gene_annot} \\
      -r {reference} \\
      -1 {fastq1} \\
      -2 {fastq2} \\
      -p {output_dir}""".format(
            genome_build=global_conf.global_get(ini_section, 'genome_build'),
            gene_annot=global_conf.global_get(ini_section, 'gene_annot'),
            reference=global_conf.global_get(ini_section, 'reference'),
            threads=global_conf.global_get(ini_section, 'threads', param_type='posint'),
            options=global_conf.global_get(ini_section, 'options'),
            fastq1=",".join(fastq1 for fastq1 in fastqs1),
            fastq2=",".join(fastq2 for fastq2 in fastqs2),
            output_dir=output_dir,
        ),
    )
