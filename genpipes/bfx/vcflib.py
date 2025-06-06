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
from ..core.job import Job

def vcfsamplediff(input_normal, input_tumor, input_vcf, output):
    return Job(
        [input_vcf],
        [output],
        [
            ['vcflib_vcfsamplediff', 'module_vcflib']
        ],
        command="""\
vcfsamplediff \\
  STATUS \\
  {input_normal} {input_tumor} \\
  {input_vcf} \\
  {output}""".format(
        input_normal=input_normal,
        input_tumor=input_tumor,
        input_vcf=input_vcf if input_vcf else " - ",
        output=" \\\n  > " + output if output else ""
        )
    )

def vcffilter(input_vcf, output_vcf, options):
    return Job(
        [input_vcf],
        [output_vcf],
        [
            ['vcflib_vcffilter', 'module_vcflib']
        ],
        command="""\
vcffilter \\
  {options} \\
  {input_vcf} \\
  {output_vcf}""".format(
        options=options,
        input_vcf=input_vcf if input_vcf else "",
        output_vcf=" \\\n  > " + output_vcf if output_vcf else ""
        )
    )

def vcfstreamsort(input_vcf, output_vcf):
    return Job(
        [input_vcf],
        [output_vcf],
        [
            ['vcflib_vcfstreamsort', 'module_vcflib']
        ],
        command="""\
vcfstreamsort \\
  {input_vcf} \\
  {output_vcf}""".format(
        input_vcf=input_vcf if input_vcf else "",
        output_vcf=" \\\n  > " + output_vcf if output_vcf else ""
        )
    )
