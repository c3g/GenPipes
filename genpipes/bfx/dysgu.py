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

def call(
    input_bam,
    output_vcf,
    region,
    ini_section="dysgu"):
    """
    Call structural variants with dysgu.

    :return: a job for dysgu structural variant calling
    """

    return Job(
        [input_bam],
        [output_vcf],
        [
            [ini_section, "module_dysgu"],
        ],
        command="""\
dysgu call {other_options} \\
    --procs {threads} \\
    --mode {mode} --overwrite \\
    -f vcf \\
    -o {output} \\
    {region} {region_setting} \\
    {genome_fasta} \\
    {tmp} \\
    {input}""".format(
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            threads=global_conf.global_get(ini_section, 'threads'),
            mode=global_conf.global_get(ini_section, 'mode'),
            output=output_vcf,
            region="--search " + {region} if region else "",
            region_setting=global_conf.global_get(ini_section, 'region_strategy') if region and global_conf.global_get(ini_section, 'region_strategy', required=False) else "",
            genome_fasta=global_conf.global_get(ini_section, 'genome_fasta'),
            tmp=global_conf.global_get(ini_section, 'tmp_dir'),
            input=input_bam
        )
    )
