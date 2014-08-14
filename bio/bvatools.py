#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def resolve_readset_coverage_bed(readset):
    coverage_target = config.param('bvatools_depth_of_coverage', 'coverage_targets', required=False)

    if coverage_target:
        if coverage_target == "auto":
            if readset.beds:
                return readset.beds[0]
            else:
                return None
        else:
            # Add filepath validation
            coverage_target = config.param('bvatools_depth_of_coverage', 'coverage_targets', type='filepath')
            return coverage_target
    else:
        return None

# If per RG != 0 is given there will be multiple outputs, so output is a prefix
# If per RG == 0 or undef, output is an actual file.
def basefreq(input, output, positions, per_rg):
    job = Job([input, positions], [output], [['bvatools_basefreq', 'module_java'], ['bvatools_basefreq', 'module_bvatools']])

    threads = config.param('bvatools_basefreq', 'threads', type='int')

    job.command = \
"""java {java_other_options} -Xmx{ram} -jar \$BVATOOLS_JAR \\
  basefreq \\
  --pos {positions} \\
  --bam {input}{per_rg}{threads} \\
  --out {output}""".format(
        java_other_options=config.param('bvatools_basefreq', 'java_other_options'),
        ram=config.param('bvatools_basefreq', 'ram'),
        positions=positions,
        per_rg=" \\\n  --per_rg " if per_rg else "",
        threads=" \\\n  --useIndex --threads " + str(threads) if threads > 1 else "",
        input=input,
        output=output
    )

    return job

def depth_of_coverage(input, output, coverage_bed, reference_genome="", other_options=""):
    job = Job([input, coverage_bed], [output], [['bvatools_depth_of_coverage', 'module_java'], ['bvatools_depth_of_coverage', 'module_bvatools']])

    if not reference_genome:
        reference_genome=config.param('bvatools_depth_of_coverage', 'genome_fasta', type='filepath')

    job.command = \
"""java {java_other_options} -Xmx{ram} -jar \$BVATOOLS_JAR \\
  depthofcoverage {other_options} \\
  --ref {reference_genome}{intervals} \\
  --bam {input} \\
  > {output}""".format(
        java_other_options=config.param('bvatools_depth_of_coverage', 'java_other_options'),
        ram=config.param('bvatools_depth_of_coverage', 'ram'),
        other_options=other_options,
        reference_genome=reference_genome,
        intervals=" \\\n  --intervals " + coverage_bed if coverage_bed else "",
        input=input,
        output=output
    )

    return job

def extract_sclip(bamFile, output_prefix, flank="200"):
    job = Job(
        [bamFile], 
        [
          output_prefix + ".sc.bam", 
          output_prefix + ".scOthers.bam", 
          output_prefix + ".scPositions.txt", 
          output_prefix + ".scSequences.txt"
        ], 
        [
          ['bvatools_ratiobaf', 'module_java'], 
          ['bvatools_ratiobaf', 'module_bvatools']
        ])

    reference_dictionary = config.param('bvatools_ratiobaf', 'genome_dictionary', type='filepath')

    job.command = \
"""java {java_other_options} -Xmx{ram} -jar \$BVATOOLS_JAR \\
  extractsclip {other_options} \\
  --bam {bamFile} \\
  --flank {flank} \\
  --minSCCount {minSCCount} \\
  --minSCLength {minSCLength} \\
  --minMappingQuality {minMappingQuality} \\
  --threads {threads} \\
  --prefix {output_prefix}""".format(
        java_other_options=config.param('bvatools_extractsclip', 'java_other_options'),
        ram=config.param('bvatools_extractsclip', 'ram'),
        other_options=config.param('bvatools_extractsclip', 'other_options', required=False),
        bamFile=bamFile,
        flank=flank,
        minSCCount=config.param('bvatools_extractsclip', 'min_sclip_count'),
        minSCLength=config.param('bvatools_extractsclip', 'kmer'),
        minMappingQuality=config.param('bvatools_extractsclip', 'min_mapping_quality'),
        threads=config.param('bvatools_extractsclip', 'threads'),
        output_prefix=output_prefix
    )

    return job
    
def groupfixmate(input, output):
    job = Job([input], [output], [['bvatools_groupfixmate', 'module_java'], ['bvatools_groupfixmate', 'module_bvatools']])

    job.command = \
"""java {java_other_options} -Xmx{ram} -jar \$BVATOOLS_JAR \\
  groupfixmate \\
  --level 1 \\
  --bam {input} \\
  --out {output}""".format(
        java_other_options=config.param('bvatools_groupfixmate', 'java_other_options'),
        ram=config.param('bvatools_groupfixmate', 'ram'),
        input=input,
        output=output
    )

    return job

def ratiobaf(basefreq, output_prefix, positions):
    job = Job([basefreq, positions], [output_prefix + ".png"], [['bvatools_ratiobaf', 'module_java'], ['bvatools_ratiobaf', 'module_bvatools']])

    reference_dictionary = config.param('bvatools_ratiobaf', 'genome_dictionary', type='filepath')

    job.command = \
"""java {java_other_options} -Xmx{ram} -jar \$BVATOOLS_JAR \\
  ratiobaf {other_options} \\
  --refdict {reference_dictionary} \\
  --snppos {positions} \\
  --basefreq {basefreq} \\
  --prefix {output_prefix}""".format(
        java_other_options=config.param('bvatools_ratiobaf', 'java_other_options'),
        ram=config.param('bvatools_ratiobaf', 'ram'),
        other_options=config.param('bvatools_ratiobaf', 'other_options', required=False),
        reference_dictionary=reference_dictionary,
        positions=positions,
        basefreq=basefreq,
        output_prefix=output_prefix
    )

    return job
