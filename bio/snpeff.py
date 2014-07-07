#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def compute_effects(input, output, split=False):
    output_stats = output + ".stats.csv"
    job = Job([input], [output, output_stats], [['compute_effects', 'moduleVersion.java'], ['compute_effects', 'moduleVersion.snpeff']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$SNPEFF_HOME/snpEff.jar eff {options} \\
  -c \$SNPEFF_HOME/snpEff.config \\
  -i vcf \\
  -o vcf \\
  -csvStats \\
  -stats {output_stats} \\
  {reference_snpeff_genome} \\
  {input}{output}""".format(
        tmp_dir=config.param('compute_effects', 'tmpDir'),
        extra_java_flags=config.param('compute_effects', 'extraJavaFlags'),
        ram=config.param('compute_effects', 'ram'),
        options=config.param('compute_effects', 'options', required=False),
        output_stats=output_stats,
        reference_snpeff_genome=config.param('compute_effects', 'referenceSnpEffGenome'),
        input=input,
        output=" \\\n  > " + output if output else ""
    )

    if split:
        split_output_stats = output + ".statsFile.txt"
        split_job = Job([output_stats], [split_output_stats], [['compute_effects', 'moduleVersion.tools']])
        split_job.command = \
"""splitSnpEffStat.awk \\
  {output_stats} \\
  {output_part} \\
  {split_output_stats}""".format(
            output_stats=output_stats,
            output_part=output + ".part",
            split_output_stats=split_output_stats
        )

        job = concat_jobs([job, split_job])

    return job

def snpsift_annotate(input, output):

    job = Job([input], [output], [['snpsift_annotate', 'moduleVersion.java'], ['snpsift_annotate', 'moduleVersion.snpeff']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$SNPEFF_HOME/SnpSift.jar annotate \\
  {db_snp} \\
  {input}{output}""".format(
        tmp_dir=config.param('snpsift_annotate', 'tmpDir'),
        extra_java_flags=config.param('snpsift_annotate', 'extraJavaFlags'),
        ram=config.param('snpsift_annotate', 'ram'),
        db_snp=config.param('snpsift_annotate', 'dbSnp', type='filepath'),
        input=input,
        output=" \\\n  > " + output if output else ""
    )

    return job

def snpsift_dbnsfp(input, output):

    job = Job([input], [output], [['snpsift_dbnsfp', 'moduleVersion.java'], ['snpsift_dbnsfp', 'moduleVersion.snpeff']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$SNPEFF_HOME/SnpSift.jar dbnsfp \\
  -v {db_nsfp} \\
  {input}{output}""".format(
        tmp_dir=config.param('snpsift_dbnsfp', 'tmpDir'),
        extra_java_flags=config.param('snpsift_dbnsfp', 'extraJavaFlags'),
        ram=config.param('snpsift_dbnsfp', 'ram'),
        db_nsfp=config.param('snpsift_dbnsfp', 'dbNSFP', type='filepath'),
        input=input,
        output=" \\\n  > " + output if output else ""
    )

    return job

