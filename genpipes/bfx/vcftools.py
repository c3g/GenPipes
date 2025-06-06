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

# Python Standard Modules

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def annotate_mappability(input, output):
    return Job(
        [input],
        [output],
        [
            ['annotate_mappability', 'module_vcftools'],
            ['annotate_mappability', 'module_htslib']
        ],
        command="""\
vcf-annotate \\
  -d key=INFO,ID=MIL,Number=1,Type=String,Description='Mappability annotation. 300IS 40SD 1SHI. HC = too high coverage (>400), LC = too low coverage (<50), MQ = too low mean mapQ (<20), ND = no data at the position' \\
  -c CHROM,FROM,TO,INFO/MIL \\
  -a {annotations} \\
  {input}{output}""".format(
        annotations=global_conf.global_get('annotate_mappability', 'genome_mappability_bed_indexed', param_type='filepath'),
        input=input,
        output=" \\\n  > " + output if output else ""
        ),
        removable_files=[output]
    )

def missing_indv(input, output):
    out = output + ".imiss"
    return Job(
        [input],
        [out],
        [
            ['vcftools_missing_indv', 'module_vcftools'],
            ['vcftools_missing_indv', 'module_htslib']
        ],
        command="""\
vcftools {options} --missing-indv \\
  --gzvcf {input} \\
  --out {output}""".format(
        options=global_conf.global_get('vcftools_missing_indv', 'options'),
        input=input,
        output=output,
        ),
        removable_files=[output]
    )

def depth(input, output):
    out = output + ".idepth"
    return Job(
        [input],
        [out],
        [
            ['vcftools_depth', 'module_vcftools'],
            ['vcftools_depth', 'module_htslib']
        ],
        command="""\
vcftools {options} --depth \\
  --gzvcf {input} \\
  --out {output}""".format(
        options=global_conf.global_get('vcftools_depth', 'options'),
        input=input,
        output=output,
        ),
        removable_files=[output]
    )
