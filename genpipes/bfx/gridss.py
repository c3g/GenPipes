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

# GenPipes Modules
from ..core.config import global_conf
from ..core.job import Job

def paired_somatic(
    input_normal,
    input_tumor,
    normal_name,
    tumor_name,
    gridss_vcf,
    gridss_bam
    ):

    return Job(
        [
            input_normal,
            input_tumor
        ],
        [
            gridss_vcf,
            gridss_bam
        ],
        [
            ['gridss_paired_somatic', 'module_java'],
            ['gridss_paired_somatic', 'module_R'],
            ['gridss_paired_somatic', 'module_samtools'],
            ['gridss_paired_somatic', 'module_bwa'],
            ['gridss_paired_somatic', 'module_gridss']
        ],
        command="""\
$GRIDSS_HOME/gridss \\
    -t {threads} \\
    --steps preprocess,assemble,call \\
    -b {blacklist_bed} \\
    -r {reference_sequence} \\
    -w {outdir} \\
    -o {gridss_vcf} \\
    -a {gridss_bam} \\
    -j $GRIDSS_JAR \\
    --jvmheap {ram} \\
    {other_options} \\
    --labels {normal},{tumor} \\
    {input_normal} \\
    {input_tumor}""".format(
            threads=global_conf.global_get('gridss_paired_somatic', 'threads', param_type='int'),
            blacklist_bed=global_conf.global_get('gridss_paired_somatic', 'blacklist_bed', param_type='filepath'),
            reference_sequence=global_conf.global_get('gridss_paired_somatic', 'genome_fasta', param_type='filepath'),
            outdir=os.path.dirname(gridss_vcf),
            gridss_vcf=gridss_vcf,
            gridss_bam=gridss_bam,
            ram=global_conf.global_get('gridss_paired_somatic', 'ram'),
            other_options=global_conf.global_get('gridss_paired_somatic', 'other_options', param_type='string'),
            normal=normal_name,
            tumor=tumor_name,
            input_normal=input_normal,
            input_tumor=input_tumor
        )
    )
