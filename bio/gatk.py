#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def base_recalibrator(input, output):

    job = Job([input], [output], [['gatk_base_recalibrator', 'module_java'], ['gatk_base_recalibrator', 'module_gatk']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar \$GATK_JAR \\
  --analysis_type BaseRecalibrator \\
  --num_cpu_threads_per_data_thread {threads} \\
  --input_file {input} \\
  --reference_sequence {reference_sequence} \\
  --knownSites {known_sites} \\
  --out {output}""".format(
        tmp_dir=config.param('gatk_base_recalibrator', 'tmp_dir'),
        java_other_options=config.param('gatk_base_recalibrator', 'java_other_options'),
        ram=config.param('gatk_base_recalibrator', 'ram'),
        threads=config.param('gatk_base_recalibrator', 'threads', type='int'),
        input=input,
        reference_sequence=config.param('gatk_base_recalibrator', 'genome_fasta', type='filepath'),
        known_sites=config.param('gatk_base_recalibrator', 'known_sites', type='filepath'),
        output=output
    )

    return job

def callable_loci(input, output, summary):

    job = Job([input], [output, summary], [['gatk_callable_loci', 'module_java'], ['gatk_callable_loci', 'module_gatk']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar \$GATK_JAR \\
  --analysis_type CallableLoci {other_options} \\
  --input_file {input} \\
  --reference_sequence {reference_sequence} \\
  --summary {summary} \\
  --out {output}""".format(
        tmp_dir=config.param('gatk_callable_loci', 'tmp_dir'),
        java_other_options=config.param('gatk_callable_loci', 'java_other_options'),
        ram=config.param('gatk_callable_loci', 'ram'),
        other_options=config.param('gatk_callable_loci', 'other_options'),
        input=input,
        reference_sequence=config.param('gatk_callable_loci', 'genome_fasta', type='filepath'),
        summary=summary,
        output=output
    )

    return job

def cat_variants(variants, output):

    job = Job(variants, [output], [['gatk_cat_variants', 'module_java'], ['gatk_cat_variants', 'module_gatk']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -cp \$GATK_JAR \\
  org.broadinstitute.sting.tools.CatVariants {options} \\
  --reference {reference}{variants} \\
  --outputFile {output}""".format(
        tmp_dir=config.param('gatk_cat_variants', 'tmp_dir'),
        java_other_options=config.param('gatk_cat_variants', 'java_other_options'),
        ram=config.param('gatk_cat_variants', 'ram'),
        options=config.param('gatk_cat_variants', 'options'),
        reference=config.param('gatk_cat_variants', 'genome_fasta', type='filepath'),
        variants="".join(" \\\n  --variant " + variant for variant in variants),
        output=output
    )

    return job

def depth_of_coverage(input, output_prefix, intervals=""):

    job = Job([input], [output_prefix + ".sample_summary"], [['gatk_depth_of_coverage', 'module_java'], ['gatk_depth_of_coverage', 'module_gatk']])

    summary_coverage_thresholds = sorted(config.param('gatk_depth_of_coverage', 'summary_coverage_thresholds', type='list'), key=int)

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar \$GATK_JAR \\
  --analysis_type DepthOfCoverage --omitDepthOutputAtEachBase --logging_level ERROR \\
  --reference_sequence {reference_sequence} \\
  --input_file {input} \\
  --out {output_prefix}{intervals}{summary_coverage_thresholds} \\
  --start 1 --stop {highest_summary_coverage_threshold} \\
  --nBins {nbins} \\
  --downsampling_type NONE""".format(
        tmp_dir=config.param('gatk_depth_of_coverage', 'tmp_dir'),
        java_other_options=config.param('gatk_depth_of_coverage', 'java_other_options'),
        ram=config.param('gatk_depth_of_coverage', 'ram'),
        reference_sequence=config.param('gatk_depth_of_coverage', 'genome_fasta', type='filepath'),
        input=input,
        output_prefix=output_prefix,
        intervals=" \\\n  --intervals " + intervals if intervals else "",
        summary_coverage_thresholds="".join(" \\\n  --summaryCoverageThreshold " + summary_coverage_threshold for summary_coverage_threshold in summary_coverage_thresholds),
        highest_summary_coverage_threshold=summary_coverage_thresholds[-1],
        nbins=int(summary_coverage_thresholds[-1]) - 1
    )

    return job

def genotype_gvcfs(variants, output):

    job = Job(variants, [output], [['gatk_genotype_gvcfs', 'module_java'], ['gatk_genotype_gvcfs', 'module_gatk']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar \$GATK_JAR \\
  --analysis_type GenotypeGVCFs {options} \\
  --reference_sequence {reference_sequence}{variants} \\
  --out {output}""".format(
        tmp_dir=config.param('gatk_genotype_gvcfs', 'tmp_dir'),
        java_other_options=config.param('gatk_genotype_gvcfs', 'java_other_options'),
        ram=config.param('gatk_genotype_gvcfs', 'ram'),
        options=config.param('gatk_genotype_gvcfs', 'options'),
        reference_sequence=config.param('gatk_genotype_gvcfs', 'genome_fasta', type='filepath'),
        variants="".join(" \\\n  --variant " + variant for variant in variants),
        output=output
    )

    return job

def haplotype_caller(input, output, intervals=[], exclude_intervals=[]):

    job = Job([input], [output], [['gatk_haplotype_caller', 'module_java'], ['gatk_haplotype_caller', 'module_gatk']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar \$GATK_JAR \\
  --analysis_type HaplotypeCaller {options} \\
  --reference_sequence {reference_sequence} \\
  --input_file {input} \\
  --out {output}{intervals}{exclude_intervals}""".format(
        tmp_dir=config.param('gatk_haplotype_caller', 'tmp_dir'),
        java_other_options=config.param('gatk_haplotype_caller', 'java_other_options'),
        ram=config.param('gatk_haplotype_caller', 'ram'),
        options=config.param('gatk_haplotype_caller', 'options'),
        reference_sequence=config.param('gatk_haplotype_caller', 'genome_fasta', type='filepath'),
        input=input,
        output=output,
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
    )

    return job

def indel_realigner(input, output, target_intervals, intervals=[], exclude_intervals=[]):

    job = Job([input], [output], [['gatk_indel_realigner', 'module_java'], ['gatk_indel_realigner', 'module_gatk']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar \$GATK_JAR \\
  --analysis_type IndelRealigner {other_options} \\
  --reference_sequence {reference_sequence} \\
  --input_file {input} \\
  --targetIntervals {target_intervals} \\
  --out {output}{intervals}{exclude_intervals} \\
  --maxReadsInMemory {max_reads_in_memory}""".format(
        tmp_dir=config.param('gatk_indel_realigner', 'tmp_dir'),
        java_other_options=config.param('gatk_indel_realigner', 'java_other_options'),
        ram=config.param('gatk_indel_realigner', 'ram'),
        other_options=config.param('gatk_indel_realigner', 'other_options'),
        reference_sequence=config.param('gatk_indel_realigner', 'genome_fasta', type='filepath'),
        input=input,
        target_intervals=target_intervals,
        output=output,
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals),
        max_reads_in_memory=config.param('gatk_indel_realigner', 'max_reads_in_memory')
    )

    return job

def print_reads(input, output, base_quality_score_recalibration):

    job = Job([input], [output], [['gatk_print_reads', 'module_java'], ['gatk_print_reads', 'module_gatk']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar \$GATK_JAR \\
  --analysis_type PrintReads \\
  --num_cpu_threads_per_data_thread {threads} \\
  --input_file {input} \\
  --reference_sequence {reference_sequence} \\
  --BQSR {base_quality_score_recalibration} \\
  --out {output}""".format(
        tmp_dir=config.param('gatk_print_reads', 'tmp_dir'),
        java_other_options=config.param('gatk_print_reads', 'java_other_options'),
        ram=config.param('gatk_print_reads', 'ram'),
        threads=config.param('gatk_print_reads', 'threads', type='int'),
        input=input,
        reference_sequence=config.param('gatk_print_reads', 'genome_fasta', type='filepath'),
        base_quality_score_recalibration=base_quality_score_recalibration,
        output=output
    )

    return job


def realigner_target_creator(input, output, intervals=[], exclude_intervals=[]):

    job = Job([input], [output], [['gatk_realigner_target_creator', 'module_java'], ['gatk_realigner_target_creator', 'module_gatk']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar \$GATK_JAR \\
  --analysis_type RealignerTargetCreator {other_options} \\
  --reference_sequence {reference_sequence} \\
  --input_file {input} \\
  --out {output}{intervals}{exclude_intervals}""".format(
        tmp_dir=config.param('gatk_realigner_target_creator', 'tmp_dir'),
        java_other_options=config.param('gatk_realigner_target_creator', 'java_other_options'),
        ram=config.param('gatk_realigner_target_creator', 'ram'),
        other_options=config.param('gatk_realigner_target_creator', 'other_options'),
        reference_sequence=config.param('gatk_realigner_target_creator', 'genome_fasta', type='filepath'),
        input=input,
        output=output,
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
    )

    return job
