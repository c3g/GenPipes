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

# GenPipes Modules
from ..core.config import global_conf
from ..core.job import Job

def genotyper(
    input_tumor,
    input_normal,
    input_vcf,
    output_vcf,
    ini_section='wham_call_sv'
    ):

    if input_tumor is not None:
        inputs=[input_tumor, input_normal]        
    else:
        inputs=[input_normal]
        
    return Job(
        inputs,
        [output_vcf],
        [
            [ini_section, 'module_python'],
        ],
        command="""\
svtyper --max_reads 5000 \\
  -T {ref_fastq} \\
  -B {input_bam} \\
  {input_vcf} \\
  {output_vcf}""".format(
        ref_fastq=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
        input_bam=",".join([input for input in inputs]),
        input_vcf="-i " + input_vcf if input_vcf else "",
        output_vcf="> " + output_vcf if output_vcf else "",
        )
    )
