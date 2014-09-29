#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def annotate_mappability(input, output):
    job = Job([input], [output], [['annotate_mappability', 'module_vcftools'], ['annotate_mappability', 'module_tabix']])

    job.command = """\
vcf-annotate \\
  -d key=INFO,ID=MIL,Number=1,Type=String,Description='Mappability annotation. 300IS 40SD 1SHI. HC = to high coverage (>400), LC = to high coverage (<50), MQ = to low mean mapQ (<20), ND = no data at the position' \\
  -c CHROM,FROM,TO,INFO/MIL \\
  -a {annotations} \\
  {input}{output}""".format(
        annotations=config.param('annotate_mappability', 'genome_mappability_bed_indexed', type='filepath'),
        input=input,
        output=" \\\n  > " + output if output else ""
    )

    return job
