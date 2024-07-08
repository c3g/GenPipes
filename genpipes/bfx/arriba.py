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

def run(fastqs1, fastqs2, output_dir):
    output_file = os.path.join(output_dir, "fusions.tsv")
    output_log = os.path.join(output_dir, "Log.final.out")
    return Job(
        fastqs1,
        [output_file, output_log],
        [
            ['run_arriba', 'module_arriba'],
            ['run_arriba', 'module_star'],
	    ['run_arriba', 'module_samtools']
        ],

        command="""\
$ARRIBA_HOME/./run_arriba.sh \\
      {genome_build} \\
      {gene_annot} \\
      {reference} \\
      {blacklist} \\
      {known_fusions} \\
      {protein_domains} \\
      {threads} \\
      {fastq1} \\
      {fastq2}""".format(
            genome_build=global_conf.global_get('run_arriba', 'genome_build'),
            gene_annot=global_conf.global_get('run_arriba', 'gene_annot'),
            reference=global_conf.global_get('run_arriba', 'reference'),
            blacklist=global_conf.global_get('run_arriba', 'blacklist'),
            known_fusions=global_conf.global_get('run_arriba', 'known_fusions'),
            protein_domains=global_conf.global_get('run_arriba', 'protein_domains'),
            threads=global_conf.global_get('run_arriba', 'threads', param_type='posint'),
            options=global_conf.global_get('run_arriba', 'options'),
            fastq1=",".join(fastq1 for fastq1 in fastqs1),
            fastq2=",".join(fastq2 for fastq2 in fastqs2),
        ),
    )
