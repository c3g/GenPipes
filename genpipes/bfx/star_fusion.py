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

def run(fastqs1, fastqs2, output_dir):
    output_file = os.path.join(output_dir, "star-fusion.fusion_predictions.abridged.coding_effect.tsv")
    jxt_file = os.path.join(output_dir, "Chimeric.out.junction")
    star_file = os.path.join(output_dir, "Aligned.out.bam")
    log_file = os.path.join(output_dir, "Log.final.out")
    
    return Job(
        fastqs1,
        [output_file,jxt_file,star_file,log_file],
        [
            ['run_star_fusion','module_perl'],
            ['run_star_fusion','module_star'],
            ['run_star_fusion','module_samtools'],
            ['run_star_fusion','module_star_fusion'],
            #['run_star_fusion', 'module_gcc']
        ],

        command="""\
    STAR-Fusion --CPU {threads} {options} \\
        --genome_lib_dir {genome_build} \\
        --left_fq {fastq1} \\
        --right_fq {fastq2} \\
        --output_dir {output_dir}""".format(
            genome_build=global_conf.global_get('run_star_fusion', 'genome_build'),
            threads=global_conf.global_get('run_star_fusion', 'threads', param_type='posint'),
            options=global_conf.global_get('run_star_fusion', 'options'),
            fastq1=",".join(fastq1 for fastq1 in fastqs1),
            fastq2=",".join(fastq2 for fastq2 in fastqs2),
            output_dir=output_dir,
        ),
    )
