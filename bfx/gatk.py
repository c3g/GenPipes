#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def base_recalibrator(input, output):

    return Job(
        [input],
        [output],
        [
            ['gatk_base_recalibrator', 'module_java'],
            ['gatk_base_recalibrator', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
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
        known_sites=config.param('gatk_base_recalibrator', 'known_variants', type='filepath'),
        output=output
        ),
        removable_files=[output]
    )

def callable_loci(input, output, summary):

    return Job(
        [input],
        [output, summary],
        [
            ['gatk_callable_loci', 'module_java'],
            ['gatk_callable_loci', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
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
    )

def cat_variants(variants, output):

    return Job(
        variants,
        [output],
        [
            ['gatk_cat_variants', 'module_java'],
            ['gatk_cat_variants', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -cp $GATK_JAR \\
  org.broadinstitute.gatk.tools.CatVariants {options} \\
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
    )


def combine_variants(variants, output):

    return Job(
        variants,
        [output],
        [
            ['gatk_combine_variants', 'module_java'],
            ['gatk_combine_variants', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type CombineVariants  \\
  --reference_sequence {reference}{variants} \\
  --out {output}""".format(
        tmp_dir=config.param('gatk_combine_variants', 'tmp_dir'),
        java_other_options=config.param('gatk_combine_variants', 'java_other_options'),
        ram=config.param('gatk_combine_variants', 'ram'),
        reference=config.param('gatk_combine_variants', 'genome_fasta', type='filepath'),
        variants="".join(" \\\n  --variant:V" + str(idx) + " " + variant for idx,variant in enumerate(variants)),
        output=output
        )
    )


def depth_of_coverage(input, output_prefix, intervals):


    summary_coverage_thresholds = sorted(config.param('gatk_depth_of_coverage', 'summary_coverage_thresholds', type='list'), key=int)

    return Job(
        [input], 
        [
          output_prefix + ".sample_summary",
          output_prefix + ".sample_cumulative_coverage_counts",
          output_prefix + ".sample_cumulative_coverage_proportions",
          output_prefix + ".sample_interval_statistics",
          output_prefix + ".sample_interval_summary",
          output_prefix + ".sample_statistics"
        ],  
        [
          ['gatk_depth_of_coverage', 'module_java'], 
          ['gatk_depth_of_coverage', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
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
    )

def genotype_gvcf(variants, output, options):

    return Job(
        variants,
        [output],
        [
            ['gatk_genotype_gvcf', 'module_java'],
            ['gatk_genotype_gvcf', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type GenotypeGVCFs {options} \\
  --disable_auto_index_creation_and_locking_when_reading_rods \\
  --reference_sequence {reference_sequence}{variants} \\
  --out {output}""".format(
        tmp_dir=config.param('gatk_genotype_gvcf', 'tmp_dir'),
        java_other_options=config.param('gatk_genotype_gvcf', 'java_other_options'),
        ram=config.param('gatk_genotype_gvcf', 'ram'),
        options=options,
        reference_sequence=config.param('gatk_genotype_gvcf', 'genome_fasta', type='filepath'),
        variants="".join(" \\\n  --variant " + variant for variant in variants),
        output=output
        )
    )

def haplotype_caller(input, output, intervals=[], exclude_intervals=[]):

    return Job(
        [input],
        [output],
        [
            ['gatk_haplotype_caller', 'module_java'],
            ['gatk_haplotype_caller', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
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
    )

def mutect(inputNormal, inputTumor, outputStats, outputVCF, intervals=[], exclude_intervals=[]):
    cosmic = config.param('gatk_mutect', 'cosmic', type='filepath', required=False)
    # if set add arg prefix
    if cosmic :
        cosmic = " --cosmic " + cosmic

    return Job(
        [inputNormal, inputTumor],
        [outputStats, outputVCF],
        [
            ['gatk_mutect', 'module_java'],
            ['gatk_mutect', 'module_mutect']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $MUTECT_JAR \\
  --analysis_type MuTect {options} \\
  --reference_sequence {reference_sequence} \\
  --dbsnp {known_sites}{cosmic} \\
  --input_file:normal {inputNormal} \\
  --input_file:tumor {inputTumor} \\
  --out {outputStats} \\
  --vcf {outputVCF}{intervals}{exclude_intervals}""".format(
        tmp_dir=config.param('gatk_mutect', 'tmp_dir'),
        java_other_options=config.param('gatk_mutect', 'java_other_options'),
        ram=config.param('gatk_mutect', 'ram'),
        options=config.param('gatk_mutect', 'options'),
        reference_sequence=config.param('gatk_mutect', 'genome_fasta', type='filepath'),
        known_sites=config.param('gatk_mutect', 'known_variants', type='filepath'),
        cosmic=cosmic,
        inputNormal=inputNormal,
        inputTumor=inputTumor,
        outputStats=outputStats,
        outputVCF=outputVCF,
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
        )
    )

def indel_realigner(input, output, target_intervals, intervals=[], exclude_intervals=[]):

    return Job(
        [input],
        [output],
        [
            ['gatk_indel_realigner', 'module_java'],
            ['gatk_indel_realigner', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
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
    )

def print_reads(input, output, base_quality_score_recalibration):

    return Job(
        [input],
        [output, re.sub("\.([sb])am$", ".\\1ai", output)],
        [
            ['gatk_print_reads', 'module_java'],
            ['gatk_print_reads', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
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
    )

def realigner_target_creator(input, output, intervals=[], exclude_intervals=[]):

    return Job(
        [input],
        [output],
        [
            ['gatk_realigner_target_creator', 'module_java'],
            ['gatk_realigner_target_creator', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
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
    )


def combine_gvcf(inputs, output, intervals=[], exclude_intervals=[]):

    return Job(
        inputs,
        [output],
        [
            ['gatk_combine_gvcf', 'module_java'],
            ['gatk_combine_gvcf', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type CombineGVCFs {other_options} \\
  --reference_sequence {reference_sequence} \\
  {input} \\
  --out {output}{intervals}{exclude_intervals}""".format(
        tmp_dir=config.param('gatk_combine_gvcf', 'tmp_dir'),
        java_other_options=config.param('gatk_combine_gvcf', 'java_other_options'),
        ram=config.param('gatk_combine_gvcf', 'ram'),
        other_options=config.param('gatk_combine_gvcf', 'other_options',required=False),
        reference_sequence=config.param('gatk_combine_gvcf', 'genome_fasta', type='filepath'),
        input="".join(" \\\n  --variant " + input for input in inputs),
        output=output,
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
        )
    )

def variant_recalibrator(variants, other_options, recal_output, tranches_output, R_output):

    return Job(
        variants,
        [recal_output, tranches_output, R_output],
        [
            ['gatk_variant_recalibrator', 'module_java'],
            ['gatk_variant_recalibrator', 'module_gatk'],
            ['gatk_variant_recalibrator', 'module_R']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type VariantRecalibrator {options} \\
  --disable_auto_index_creation_and_locking_when_reading_rods \\
  --reference_sequence {reference_sequence}{variants} \\
  {other_options} \\
  --recal_file {recal_output} \\
  --tranches_file {tranches_output} \\
  --rscript_file {R_output}""".format(
        tmp_dir=config.param('gatk_variant_recalibrator', 'tmp_dir'),
        java_other_options=config.param('gatk_variant_recalibrator', 'java_other_options'),
        ram=config.param('gatk_variant_recalibrator', 'ram'),
        options=config.param('gatk_variant_recalibrator', 'options'),
        reference_sequence=config.param('gatk_variant_recalibrator', 'genome_fasta', type='filepath'),
        variants="".join(" \\\n  -input " + variant for variant in variants),
        other_options=other_options,
        recal_output=recal_output,
        tranches_output=tranches_output,
        R_output=R_output
        ),
        removable_files=[recal_output, tranches_output, R_output]
    )
        
def apply_recalibration(variants, recal_input, tranches_input, other_options, apply_recal_output):

    return Job(
        [variants, recal_input, tranches_input],
        [apply_recal_output],
        [
            ['gatk_apply_recalibration', 'module_java'],
            ['gatk_apply_recalibration', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type ApplyRecalibration {options} \\
  --disable_auto_index_creation_and_locking_when_reading_rods \\
  --reference_sequence {reference_sequence} \\
  -input {variants} \\
  {other_options} \\
  --tranches_file {tranches_input} \\
  --recal_file {recal_input} \\
  --out {output}""".format(
        tmp_dir=config.param('gatk_apply_recalibration', 'tmp_dir'),
        java_other_options=config.param('gatk_apply_recalibration', 'java_other_options'),
        ram=config.param('gatk_apply_recalibration', 'ram'),
        options=config.param('gatk_apply_recalibration', 'options'),
        reference_sequence=config.param('gatk_apply_recalibration', 'genome_fasta', type='filepath'),
        variants=variants,
        other_options=other_options,
        recal_input=recal_input,
        tranches_input=tranches_input,
        output=apply_recal_output
        )
    )

