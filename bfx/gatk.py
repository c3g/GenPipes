################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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
import os

# MUGQIC Modules
import core
from core.job import Job
from core.config import global_config_parser
from . import gatk4
from . import picard
from . import picard2

config = core.config.global_config_parser

def base_recalibrator(input, output, intervals):
    if intervals:
        inputs = [input, intervals]
    
    else:
        inputs = [input]
        
    if global_config_parser.param('gatk_base_recalibrator', 'module_gatk').split("/")[2] >= "4":
        return gatk4.base_recalibrator(input, output, intervals)
    else:
        return Job(
            inputs,
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
        tmp_dir=global_config_parser.param('gatk_base_recalibrator', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_base_recalibrator', 'java_other_options'),
        options=global_config_parser.param('gatk_base_recalibrator', 'options'),
        ram=global_config_parser.param('gatk_base_recalibrator', 'ram'),
        threads=global_config_parser.param('gatk_base_recalibrator', 'threads', param_type='int'),
        input=input,
        intervals=" \\\n  --intervals " + intervals if intervals else "",
        reference_sequence=global_config_parser.param('gatk_base_recalibrator', 'genome_fasta', param_type='filepath'),
        known_dbsnp=global_config_parser.param('gatk_base_recalibrator', 'known_dbsnp', param_type='filepath'),
        known_gnomad=global_config_parser.param('gatk_base_recalibrator', 'known_gnomad', param_type='filepath'),
        known_mills=global_config_parser.param('gatk_base_recalibrator', 'known_mills', param_type='filepath'),
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
        tmp_dir=global_config_parser.param('gatk_callable_loci', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_callable_loci', 'java_other_options'),
        ram=global_config_parser.param('gatk_callable_loci', 'ram'),
        other_options=global_config_parser.param('gatk_callable_loci', 'other_options'),
        input=input,
        reference_sequence=global_config_parser.param('gatk_callable_loci', 'genome_fasta', param_type='filepath'),
        summary=summary,
        output=output
        )
    )

def cat_variants(variants, output=None):
    if global_config_parser.param('gatk_cat_variants', 'module_gatk').split("/")[2] >= "4":
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
        tmp_dir=global_config_parser.param('gatk_cat_variants', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_cat_variants', 'java_other_options'),
        ram=global_config_parser.param('gatk_cat_variants', 'ram'),
        options=global_config_parser.param('gatk_cat_variants', 'options'),
        reference=global_config_parser.param('gatk_cat_variants', 'genome_fasta', param_type='filepath'),
        variants="".join(" \\\n  --variant " + variant for variant in variants),
        output=output
        )
    )


def combine_variants(variants, output):
    if global_config_parser.param('gatk_combine_variants', 'module_gatk').split("/")[2] >= "4":
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
        tmp_dir=global_config_parser.param('gatk_combine_variants', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_combine_variants', 'java_other_options'),
        ram=global_config_parser.param('gatk_combine_variants', 'ram'),
        reference=global_config_parser.param('gatk_combine_variants', 'genome_fasta', param_type='filepath'),
        variants="".join(" \\\n  --variant:V" + str(idx) + " " + variant for idx,variant in enumerate(variants)),
        output=output
        )
    )


def depth_of_coverage(input, output_prefix, intervals):

    summary_coverage_thresholds = sorted(global_config_parser.param('gatk_depth_of_coverage', 'summary_coverage_thresholds', param_type='list'), key=int)

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
        tmp_dir=global_config_parser.param('gatk_depth_of_coverage', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_depth_of_coverage', 'java_other_options'),
        ram=global_config_parser.param('gatk_depth_of_coverage', 'ram'),
        reference_sequence=global_config_parser.param('gatk_depth_of_coverage', 'genome_fasta', param_type='filepath'),
        input=input,
        output_prefix=output_prefix,
        intervals=" \\\n  --intervals " + intervals if intervals else "",
        summary_coverage_thresholds="".join(" \\\n  --summaryCoverageThreshold " + summary_coverage_threshold for summary_coverage_threshold in summary_coverage_thresholds),
        highest_summary_coverage_threshold=summary_coverage_thresholds[-1],
        nbins=int(summary_coverage_thresholds[-1]) - 1
        )
    )

def genotype_gvcf(variants, output, options):
    if global_config_parser.param('gatk_genotype_gvcf', 'module_gatk').split("/")[2] >= "4":
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
        tmp_dir=global_config_parser.param('gatk_genotype_gvcf', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_genotype_gvcf', 'java_other_options'),
        ram=global_config_parser.param('gatk_genotype_gvcf', 'ram'),
        options=options,
        reference_sequence=global_config_parser.param('gatk_genotype_gvcf', 'genome_fasta', param_type='filepath'),
        variants="".join(" \\\n  --variant " + variant for variant in variants),
        output=output
        )
    )

def haplotype_caller(
    inputs,
    output,
    intervals=[],
    exclude_intervals=[],
    interval_list=None
    ):
    interval_padding = config.param('gatk_haplotype_caller', 'interval_padding')
    if not isinstance(inputs, list):
        inputs = [inputs]

    inputs_list = inputs.copy()
    if not interval_list is None:
       inputs_list.extend([interval_list])

    if global_config_parser.param('gatk_haplotype_caller', 'module_gatk').split("/")[2] >= "4":
        return gatk4.haplotype_caller(
            inputs,
            output,
            intervals=intervals,
            exclude_intervals=exclude_intervals,
            interval_list=interval_list
        )
    else:
        return Job(
            inputs_list,
            [output, output + ".tbi"],
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
  --out {output}{interval_padding} {interval_list}{intervals}{exclude_intervals}""".format(
        tmp_dir=global_config_parser.param('gatk_haplotype_caller', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_haplotype_caller', 'java_other_options'),
        ram=global_config_parser.param('gatk_haplotype_caller', 'ram'),
        options=global_config_parser.param('gatk_haplotype_caller', 'options'),
        reference_sequence=global_config_parser.param('gatk_haplotype_caller', 'genome_fasta', param_type='filepath'),
        interval_list=" --intervals " + interval_list if interval_list else "",
        interval_padding=" \\\n --interval-padding " + str(interval_padding) if interval_padding else "",
        input=" \\\n  ".join(input for input in inputs),
        output=output,
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
        )
    )

def mutect(inputNormal, inputTumor, outputStats, outputVCF, intervals=[], exclude_intervals=[]):
    cosmic = global_config_parser.param('gatk_mutect', 'cosmic', param_type='filepath', required=False)
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
        tmp_dir=global_config_parser.param('gatk_mutect', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_mutect', 'java_other_options'),
        ram=global_config_parser.param('gatk_mutect', 'ram'),
        options=global_config_parser.param('gatk_mutect', 'options'),
        reference_sequence=global_config_parser.param('gatk_mutect', 'genome_fasta', param_type='filepath'),
        known_sites=global_config_parser.param('gatk_mutect', 'known_variants', param_type='filepath'),
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
    cosmic = global_config_parser.param('gatk_mutect2', 'cosmic', param_type='filepath', required=False)
    # if set add arg prefix
    if cosmic and os.path.isfile(cosmic):
        cosmic = " --cosmic " + cosmic
    
    if interval_list:
        inputs = [inputNormal, inputTumor, interval_list]
        
    else:
        inputs = [inputNormal, inputTumor]
    
    if global_config_parser.param('gatk_mutect2', 'module_gatk').split("/")[2] >= "4":
        return gatk4.mutect2(inputNormal, normal_name, inputTumor, tumor_name, outputVCF, intervals, exclude_intervals, interval_list)
    else:
        return Job(
            inputs,
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
  --out {outputVCF}{intervals}{exclude_intervals}{interval_list}""".format(
        tmp_dir=global_config_parser.param('gatk_mutect', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_mutect2', 'java_other_options'),
        ram=global_config_parser.param('gatk_mutect2', 'ram'),
        options=global_config_parser.param('gatk_mutect2', 'options'),
        reference_sequence=global_config_parser.param('gatk_mutect2', 'genome_fasta', param_type='filepath'),
        known_sites=global_config_parser.param('gatk_mutect2', 'dbsnp', param_type='filepath'),
        cosmic=cosmic,
        inputNormal=inputNormal,
        inputTumor=inputTumor,
        outputVCF=outputVCF,
        #interval_list=" \\\n  --interval_padding 100 --intervals " + interval_list if interval_list else "",
        interval_list=" \\\n  --intervals " + interval_list if interval_list else "",
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
        )
    )

def indel_realigner(input,
                    target_intervals,
                    input2=[],
                    output=[],
                    output_dir=[],
                    output_norm_dep=[],
                    output_tum_dep=[],
                    intervals=[],
                    exclude_intervals=[],
                    optional=[]
                    ):
    
    output_dep = output_norm_dep + output_tum_dep
    output_dep.append(output)
    
    return Job(
        [input, target_intervals, input2],
        output_dep,
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
  {known_mills}{output}{intervals}{exclude_intervals} \\
  --maxReadsInMemory {max_reads_in_memory}""".format(
        tmp_dir=global_config_parser.param('gatk_indel_realigner', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_indel_realigner', 'java_other_options'),
        ram=global_config_parser.param('gatk_indel_realigner', 'ram'),
        other_options=global_config_parser.param('gatk_indel_realigner', 'other_options'),
        reference_sequence=global_config_parser.param('gatk_indel_realigner', 'genome_fasta', param_type='filepath'),
        optional="--nWayOut " + optional if optional else "",
        input=os.path.join(output_dir, input),
        input2="--input_file " + os.path.join(output_dir, input2) if input2 else "",
        target_intervals=os.path.join(output_dir, target_intervals),
        known_mills=" \\\n  --knownAlleles " + global_config_parser.param('gatk_realigner_target_creator', 'known_mills', param_type='filepath') if global_config_parser.param('gatk_realigner_target_creator', 'known_mills', param_type='filepath') else "",
        output=" \\\n  --out " + os.path.join(output_dir,output) if output else "",
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals),
        max_reads_in_memory=global_config_parser.param('gatk_indel_realigner', 'max_reads_in_memory')
        )
    )

def print_reads(input, output, base_quality_score_recalibration):
    if global_config_parser.param('gatk_print_reads', 'module_gatk').split("/")[2] >= "4":
        return gatk4.print_reads(input, output, base_quality_score_recalibration)
    else:
        return Job(
            [input, base_quality_score_recalibration],
            [output, re.sub(".bam", ".bai", output)],
            [
                ['gatk_print_reads', 'module_java'],
                ['gatk_print_reads', 'module_gatk']
            ],
        command="""\
rm -rf {output}* && \\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type PrintReads --generate_md5 \\
  -nt 1 --num_cpu_threads_per_data_thread {threads} \\
  --input_file {input} \\
  --reference_sequence {reference_sequence} \\
  --BQSR {base_quality_score_recalibration} \\
  --out {output}""".format(
        tmp_dir=global_config_parser.param('gatk_print_reads', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_print_reads', 'java_other_options'),
        ram=global_config_parser.param('gatk_print_reads', 'ram'),
        threads=global_config_parser.param('gatk_print_reads', 'threads', param_type='int'),
        input=input,
        reference_sequence=global_config_parser.param('gatk_print_reads', 'genome_fasta', param_type='filepath'),
        base_quality_score_recalibration=base_quality_score_recalibration,
        output=output
        )
    )

def realigner_target_creator(input,
                             output,
                             output_dir=[],
                             input2=[],
                             intervals=[],
                             exclude_intervals=[]
                             ):

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
  --out {output}{known_mills}{intervals}{exclude_intervals}""".format(
        tmp_dir=global_config_parser.param('gatk_realigner_target_creator', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_realigner_target_creator', 'java_other_options'),
        ram=global_config_parser.param('gatk_realigner_target_creator', 'ram'),
        other_options=global_config_parser.param('gatk_realigner_target_creator', 'other_options'),
        reference_sequence=global_config_parser.param('gatk_realigner_target_creator', 'genome_fasta', param_type='filepath'),
        input= os.path.join(output_dir, input),
        input2="--input_file " + os.path.join(output_dir, input2) if input2 else "",
        known_mills=" \\\n  --known " + global_config_parser.param('gatk_realigner_target_creator', 'known_mills', param_type='filepath') if global_config_parser.param('gatk_realigner_target_creator', 'known_mills', param_type='filepath') else "",
        output=os.path.join(output_dir, output),
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
        )
    )


def combine_gvcf(inputs, output, intervals=[], exclude_intervals=[]):

    if not isinstance(inputs, list):
        inputs=[inputs]
    
    if global_config_parser.param('gatk_combine_gvcf', 'module_gatk').split("/")[2] >= "4":
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
        tmp_dir=global_config_parser.param('gatk_combine_gvcf', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_combine_gvcf', 'java_other_options'),
        ram=global_config_parser.param('gatk_combine_gvcf', 'ram'),
        other_options=global_config_parser.param('gatk_combine_gvcf', 'other_options', required=False),
        reference_sequence=global_config_parser.param('gatk_combine_gvcf', 'genome_fasta', param_type='filepath'),
        input="".join(" \\\n  --variant " + input for input in inputs),
        output=output,
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
        )
    )

def variant_annotator(input_normal, input_tumor, input_variants, output, other_options, intervals=[], exclude_intervals=[]):

    return Job(
        [input_normal, input_tumor, input_variants],
        [output, output + ".tbi"],
        [
            ['gatk_variant_annotator', 'module_java'],
            ['gatk_variant_annotator', 'module_gatk']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GATK_JAR \\
  --analysis_type VariantAnnotator {other_options} \\
  --disable_auto_index_creation_and_locking_when_reading_rods \\
  --reference_sequence {reference_sequence} \\
  --input_file {input_normal} --input_file {input_tumor} \\
  --variant {input_variants} \\
  --out {output}{intervals}{exclude_intervals}""".format(
        tmp_dir=global_config_parser.param('gatk_variant_annotator', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_variant_annotator', 'java_other_options'),
        ram=global_config_parser.param('gatk_variant_annotator', 'ram'),
        other_options=other_options,
        reference_sequence=global_config_parser.param('gatk_variant_annotator', 'genome_fasta', param_type='filepath'),
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
    if global_config_parser.param('gatk_print_reads', 'module_gatk').split("/")[2] >= "4":
        return gatk4.combine_gvcf(variants, other_options, recal_output, tranches_output, R_output)
    else:

        if small_sample_check:
            try:
                small_sample_option = global_config_parser.param('gatk_variant_recalibrator', 'small_sample_option')
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
        tmp_dir=global_config_parser.param('gatk_variant_recalibrator', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_variant_recalibrator', 'java_other_options'),
        ram=global_config_parser.param('gatk_variant_recalibrator', 'ram'),
        options=global_config_parser.param('gatk_variant_recalibrator', 'options'),
        reference_sequence=global_config_parser.param('gatk_variant_recalibrator', 'genome_fasta', param_type='filepath'),
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
        tmp_dir=global_config_parser.param('gatk_apply_recalibration', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_apply_recalibration', 'java_other_options'),
        ram=global_config_parser.param('gatk_apply_recalibration', 'ram'),
        options=global_config_parser.param('gatk_apply_recalibration', 'options'),
        reference_sequence=global_config_parser.param('gatk_apply_recalibration', 'genome_fasta', param_type='filepath'),
        variants=variants,
        other_options=other_options,
        recal_input=recal_input,
        tranches_input=tranches_input,
        output=apply_recal_output
        )
    )


def split_n_cigar_reads(input, output, intervals=[], exclude_intervals=[], interval_list=None):
    if interval_list:
        inputs = [input]
    
    else:
        inputs = [input]
    return Job(
        inputs,
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
            tmp_dir=global_config_parser.param('gatk_split_N_trim', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_split_N_trim', 'java_other_options'),
            ram=global_config_parser.param('gatk_split_N_trim', 'ram'),
            other_options=global_config_parser.param('gatk_split_N_trim', 'other_options', required=False),
            reference_sequence=global_config_parser.param('gatk_split_N_trim', 'reference', param_type='filepath'),
            input=" \\\n  ".join(input for input in inputs),
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
            tmp_dir=global_config_parser.param('gatk_variant_filtration', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_variant_filtration', 'java_other_options'),
            ram=global_config_parser.param('gatk_variant_filtration', 'ram'),
            reference_sequence=global_config_parser.param('gatk_variant_filtration', 'genome_fasta', param_type='filepath'),
            variants=input,
            other_options=other_options,
            output=output
        )
    )

def bed2interval_list(
    dictionary,
    bed,
    output
    ):
    if global_config_parser.param('picard_bed2interval_list', 'module_gatk').split("/")[2] >= "4":
        return gatk4.bed2interval_list(
            dictionary,
            bed,
            output
        )
    elif global_config_parser.param('picard_bed2interval_list', 'module_picard').split("/")[2] < "2":
        return picard.bed2interval_list(
            dictionary,
            bed,
            output
        )
    elif global_config_parser.param('picard_bed2interval_list', 'module_picard').split("/")[2] >= "2":
        return picard2.bed2interval_list(
            dictionary,
            bed,
            output
        )
