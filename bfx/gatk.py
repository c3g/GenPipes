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
import re

# MUGQIC Modules
from core.job import *
import core.config
import gatk4

config = core.config.config

def base_recalibrator(input, output, intervals):
    if config.param('gatk_base_recalibrator', 'module_gatk').split("/")[2] >= "4":
        return gatk4.base_recalibrator(input, output, intervals)
    else:
        return Job(
            [input, intervals],
            [output],
            [
                ['gatk_base_recalibrator', 'module_java'],
                ['gatk_base_recalibrator', 'module_gatk']
            ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type BaseRecalibrator {options} \\
  -nt 1 --num_cpu_threads_per_data_thread {threads} \\
  --input_file {input} \\
  --reference_sequence {reference_sequence} {intervals} \\
  --knownSites {known_dbsnp} \\
  --knownSites {known_gnomad} \\
  --knownSites {known_mills} \\
  --out {output}""".format(
        tmp_dir=config.param('gatk_base_recalibrator', 'tmp_dir'),
        java_other_options=config.param('gatk_base_recalibrator', 'java_other_options'),
        options=config.param('gatk_base_recalibrator', 'options'),
        ram=config.param('gatk_base_recalibrator', 'ram'),
        threads=config.param('gatk_base_recalibrator', 'threads', type='int'),
        input=input,
        intervals=" \\\n  --intervals " + intervals if intervals else "",
        reference_sequence=config.param('gatk_base_recalibrator', 'genome_fasta', type='filepath'),
        known_dbsnp=config.param('gatk_base_recalibrator', 'known_dbsnp', type='filepath'),
        known_gnomad=config.param('gatk_base_recalibrator', 'known_gnomad', type='filepath'),
        known_mills=config.param('gatk_base_recalibrator', 'known_mills', type='filepath'),
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

def cat_variants(variants, output=None):
    if config.param('gatk_cat_variants', 'module_gatk').split("/")[2] >= "4":
        return gatk4.cat_variants(variants, output)
    else:
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
    if config.param('gatk_combine_variants', 'module_gatk').split("/")[2] >= "4":
        return gatk4.combine_variants(variants, output)
    else:
        return Job(
            variants,
            [output],
            [
                ['gatk_combine_variants', 'module_java'],
                ['gatk_combine_variants', 'module_gatk']
            ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type CombineVariants --genotypemergeoption UNSORTED \\
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
    if config.param('gatk_genotype_gvcf', 'module_gatk').split("/")[2] >= "4":
        return gatk4.genotype_gvcf(variants, output, options)
    else:
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

def haplotype_caller(inputs, output, intervals=[], exclude_intervals=[], interval_list=None):
    if not isinstance(inputs, list):
        inputs=[inputs, interval_list]
        
    if config.param('gatk_haplotype_caller', 'module_gatk').split("/")[2] >= "4":
        return gatk4.haplotype_caller(inputs, output, intervals, exclude_intervals, interval_list)
    else:
        return Job(
            inputs,
            [output,output+".tbi"],
            [
                ['gatk_haplotype_caller', 'module_java'],
                ['gatk_haplotype_caller', 'module_gatk']
            ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type HaplotypeCaller {options} \\
  --disable_auto_index_creation_and_locking_when_reading_rods \\
  --reference_sequence {reference_sequence} \\
  --input_file {input} \\
  --out {output}{interval_list}{intervals}{exclude_intervals}""".format(
        tmp_dir=config.param('gatk_haplotype_caller', 'tmp_dir'),
        java_other_options=config.param('gatk_haplotype_caller', 'java_other_options'),
        ram=config.param('gatk_haplotype_caller', 'ram'),
        options=config.param('gatk_haplotype_caller', 'options'),
        reference_sequence=config.param('gatk_haplotype_caller', 'genome_fasta', type='filepath'),
        interval_list=" \\\n  --interval_padding 100 --intervals " + interval_list if interval_list else "",
        input=" \\\n  ".join(input for input in inputs),
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

def mutect2(inputNormal, normal_name, inputTumor, tumor_name, outputVCF, intervals=[], exclude_intervals=[], interval_list=None):
    cosmic = config.param('gatk_mutect2', 'cosmic', type='filepath', required=False)
    # if set add arg prefix
    if cosmic :
        cosmic = " --cosmic " + cosmic
        
    if config.param('gatk_mutect2', 'module_gatk').split("/")[2] >= "4":
        return gatk4.mutect2(inputNormal, normal_name, inputTumor, tumor_name, outputVCF, intervals, exclude_intervals, interval_list)
    else:
        return Job(
            [inputNormal, inputTumor, interval_list],
            [outputVCF],
            [
                ['gatk_mutect2', 'module_java'],
                ['gatk_mutect2', 'module_gatk']
            ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type MuTect2 {options} \\
  --disable_auto_index_creation_and_locking_when_reading_rods \\
  --reference_sequence {reference_sequence} \\
  --dbsnp {known_sites}{cosmic} \\
  --input_file:normal {inputNormal} \\
  --input_file:tumor {inputTumor} \\
  --out {outputVCF}{interval_list}{intervals}{exclude_intervals}""".format(
        tmp_dir=config.param('gatk_mutect', 'tmp_dir'),
        java_other_options=config.param('gatk_mutect2', 'java_other_options'),
        ram=config.param('gatk_mutect2', 'ram'),
        options=config.param('gatk_mutect2', 'options'),
        reference_sequence=config.param('gatk_mutect2', 'genome_fasta', type='filepath'),
        known_sites=config.param('gatk_mutect2', 'dbsnp', type='filepath'),
        cosmic=cosmic,
        inputNormal=inputNormal,
        inputTumor=inputTumor,
        outputVCF=outputVCF,
        interval_list=" \\\n  --interval-padding 100 --intervals " + interval_list if interval_list else "",
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
        )
    )

def indel_realigner(input, target_intervals, input2=[], output=[], output_norm_dep=[], output_tum_dep=[], intervals=[], exclude_intervals=[], optional=[]):

    return Job(
        [input, target_intervals, input2],
        [output, output_norm_dep, output_tum_dep],
        [
            ['gatk_indel_realigner', 'module_java'],
            ['gatk_indel_realigner', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type IndelRealigner {other_options} \\
  --reference_sequence {reference_sequence} \\
  {optional} \\
  --input_file {input} \\
  {input2} \\
  --targetIntervals {target_intervals} \\
  --knownAlleles {known_mills} \\
  {output}{intervals}{exclude_intervals} \\
  --maxReadsInMemory {max_reads_in_memory}""".format(
        tmp_dir=config.param('gatk_indel_realigner', 'tmp_dir'),
        java_other_options=config.param('gatk_indel_realigner', 'java_other_options'),
        ram=config.param('gatk_indel_realigner', 'ram'),
        other_options=config.param('gatk_indel_realigner', 'other_options'),
        reference_sequence=config.param('gatk_indel_realigner', 'genome_fasta', type='filepath'),
        optional="--nWayOut " + optional if optional else "",
        input=input,
        input2="--input_file " + input2 if input2 else "",
        target_intervals=target_intervals,
        #known_indel_sites=config.param('gatk_realigner_target_creator', 'known_indel_sites', type='filepath'),
        known_mills=config.param('gatk_realigner_target_creator', 'known_mills', type='filepath'),
        #known_1000G=config.param('gatk_realigner_target_creator', 'known_1000G', type='filepath'),
        output=" \\\n  --out " + output if output else "",
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals),
        max_reads_in_memory=config.param('gatk_indel_realigner', 'max_reads_in_memory')
        )
    )

def print_reads(input, output, base_quality_score_recalibration):
    if config.param('gatk_print_reads', 'module_gatk').split("/")[2] >= "4":
        return gatk4.print_reads(input, output, base_quality_score_recalibration)
    else:
        return Job(
            [input, base_quality_score_recalibration],
            [output, re.sub("\.([sb])am$", ".\\1ai", output)],
            [
                ['gatk_print_reads', 'module_java'],
                ['gatk_print_reads', 'module_gatk']
            ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type PrintReads --generate_md5 \\
  -nt 1 --num_cpu_threads_per_data_thread {threads} \\
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

def realigner_target_creator(input, output, input2=[], intervals=[], exclude_intervals=[]):

    return Job(
        [input, input2],
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
  {input2} \\
  --known {known_mills} \\
  --out {output}{intervals}{exclude_intervals}""".format(
        tmp_dir=config.param('gatk_realigner_target_creator', 'tmp_dir'),
        java_other_options=config.param('gatk_realigner_target_creator', 'java_other_options'),
        ram=config.param('gatk_realigner_target_creator', 'ram'),
        other_options=config.param('gatk_realigner_target_creator', 'other_options'),
        reference_sequence=config.param('gatk_realigner_target_creator', 'genome_fasta', type='filepath'),
        input=input,
        input2="--input_file " + input2 if input2 else "",
        known_mills=config.param('gatk_realigner_target_creator', 'known_mills', type='filepath'),
        #known_1000G=config.param('gatk_realigner_target_creator', 'known_1000G', type='filepath'),
        output=output,
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
        )
    )


def combine_gvcf(inputs, output, intervals=[], exclude_intervals=[]):

    if not isinstance(inputs, list):
        inputs=[inputs]
    
    if config.param('gatk_combine_gvcf', 'module_gatk').split("/")[2] >= "4":
        return gatk4.combine_gvcf(inputs, output, intervals, exclude_intervals)
    else:
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
  --disable_auto_index_creation_and_locking_when_reading_rods \\
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

def variant_annotator(input_normal, input_tumor, input_variants, output, intervals=[], exclude_intervals=[]):

    return Job(
        [input_normal, input_tumor, input_variants],
        [output],
        [
            ['gatk_variant_annotator', 'module_java'],
            ['gatk_variant_annotator', 'module_gatk']
        ],
        command="""\
java {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type VariantAnnotator {other_options} \\
  --disable_auto_index_creation_and_locking_when_reading_rods \\
  --reference_sequence {reference_sequence} \\
  --input_file {input_normal} --input_file {input_tumor} \\
  --variant {input_variants} \\
  --out {output}{intervals}{exclude_intervals}""".format(
        tmp_dir=config.param('gatk_variant_annotator', 'tmp_dir'),
        java_other_options=config.param('gatk_variant_annotator', 'java_other_options'),
        ram=config.param('gatk_variant_annotator', 'ram'),
        other_options=config.param('gatk_variant_annotator', 'other_options',required=False),
        reference_sequence=config.param('gatk_variant_annotator', 'genome_fasta', type='filepath'),
        input_normal=input_normal,
        input_tumor=input_tumor,
        input_variants=input_variants,
        output=output,
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
        )
    )


def variant_recalibrator(variants, other_options, recal_output,
                         tranches_output, R_output, small_sample_check=False):
    if config.param('gatk_print_reads', 'module_gatk').split("/")[2] >= "4":
        return gatk4.combine_gvcf(variants, other_options, recal_output, tranches_output, R_output)
    else:

        if small_sample_check:
            try:
                small_sample_option = config.param('gatk_variant_recalibrator', 'small_sample_option')
            except core.config.Error:
                small_sample_option = ''
        else:
            small_sample_option = ''

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
  --tranches_file {tranches_output} {small_sample_option} \\
  --rscript_file {R_output}""".format(
        tmp_dir=config.param('gatk_variant_recalibrator', 'tmp_dir'),
        java_other_options=config.param('gatk_variant_recalibrator', 'java_other_options'),
        ram=config.param('gatk_variant_recalibrator', 'ram'),
        options=config.param('gatk_variant_recalibrator', 'options'),
        reference_sequence=config.param('gatk_variant_recalibrator', 'genome_fasta', type='filepath'),
        variants="".join(" \\\n  -input " + variant for variant in variants),
        other_options=other_options,
        #tmp_dir="--TMP_DIR " + tmp_dir if tmp_dir else "",
        recal_output=recal_output,
        tranches_output=tranches_output,
        small_sample_option=small_sample_option,
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


def split_n_cigar_reads(input, output, intervals=[], exclude_intervals=[], interval_list=None):
    return Job(
        [input],
        [output],
        [
            ['gatk_split_N_trim', 'module_java'],
            ['gatk_split_N_trim', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type SplitNCigarReads {other_options} \\
  --reference_sequence {reference_sequence} \\
  --input_file {input} \\
  --out {output}{intervals}{exclude_intervals}""".format(
            tmp_dir=config.param('gatk_split_N_trim', 'tmp_dir'),
            java_other_options=config.param('gatk_split_N_trim', 'java_other_options'),
            ram=config.param('gatk_split_N_trim', 'ram'),
            other_options=config.param('gatk_split_N_trim', 'other_options', required=False),
            reference_sequence=config.param('gatk_split_N_trim', 'reference', type='filepath'),
            input=input,
            output=output,
            intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
            interval_list=" \\\n --interval-padding 100 --intervals " + interval_list if interval_list else "",
            exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
        )
    )


def variant_filtration(input, output, other_options):
    return Job(
        [input],
        [output],
        [
            ['gatk_variant_filtration', 'module_java'],
            ['gatk_variant_filtration', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type VariantFiltration \\
  --disable_auto_index_creation_and_locking_when_reading_rods \\
  --reference_sequence {reference_sequence} \\
  --variant {variants} \\
  {other_options} \\
  --out {output}""".format(
            tmp_dir=config.param('gatk_variant_filtration', 'tmp_dir'),
            java_other_options=config.param('gatk_variant_filtration', 'java_other_options'),
            ram=config.param('gatk_variant_filtration', 'ram'),
            reference_sequence=config.param('gatk_variant_filtration', 'genome_fasta', type='filepath'),
            variants=input,
            other_options=other_options,
            output=output
        )
    )
