#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def base_recalibrator(input, output):

    job = Job([input], [output], [['baseRecalibrator', 'moduleVersion.java'], ['baseRecalibrator', 'moduleVersion.gatk']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$GATK._JAR \\
  -T {tmp_dir} \\
  -nct {threads} \\
  -I {input} \\
  -R {reference_fasta} \\
  -knownSites {known_sites} \\
  -o {output}""".format(
        tmp_dir=config.param('baseRecalibrator', 'tmpDir'),
        extra_java_flags=config.param('baseRecalibrator', 'extraJavaFlags'),
        ram=config.param('baseRecalibrator', 'ram'),
        threads=config.param('baseRecalibrator', 'threads', type='int'),
        input=input,
        reference_fasta=config.param('baseRecalibrator', 'referenceFasta', type='filepath'),
        reference_fasta=config.param('baseRecalibrator', 'knownSites', type='filepath'),
        output=output
    )

    return job
