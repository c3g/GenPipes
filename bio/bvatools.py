#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def resolve_readset_coverage_bed(readset):
    coverage_target = config.param('depth_of_coverage', 'coverageTargets', required=False)

    if coverage_target:
        if coverage_target == "auto":
            if readset.beds:
                return readset.beds[0]
            else:
                return None
        else:
            return coverage_target
    else:
        return None

def depth_of_coverage(input, output, coverage_bed, reference_genome=""):
    job = Job([input, coverage_bed], [output], [['depth_of_coverage', 'moduleVersion.java'], ['depth_of_coverage', 'moduleVersion.bvatools']])

    if not reference_genome:
        reference_genome=config.param('depth_of_coverage', 'referenceFasta', type='filepath')

    job.command = \
"""java {extra_java_flags} -Xmx{ram} -jar \$BVATOOLS_JAR \\
  depthofcoverage {extra_flags} \\
  --ref {reference_genome}{intervals} \\
  --bam {input} \\
  > {output}""".format(
        extra_java_flags=config.param('depth_of_coverage', 'extraJavaFlags'),
        ram=config.param('depth_of_coverage', 'ram'),
        extra_flags=config.param('depth_of_coverage', 'extraFlags', required=False),
        reference_genome=reference_genome,
        intervals=" \\\n  --intervals " + coverage_bed if coverage_bed else "",
        input=input,
        output=output
    )

    return job

def fix_mate_by_coordinate(input, output):
    job = Job([input], [output], [['fix_mate_by_coordinate', 'moduleVersion.java'], ['fix_mate_by_coordinate', 'moduleVersion.bvatools']])

    job.command = \
"""java {extra_java_flags} -Xmx{ram} -jar \$BVATOOLS_JAR \\
  groupfixmate \\
  --level 1 \\
  --bam {input} \\
  --out {output}""".format(
        extra_java_flags=config.param('fix_mate_by_coordinate', 'extraJavaFlags'),
        ram=config.param('fix_mate_by_coordinate', 'ram'),
        input=input,
        output=output
    )

    return job
