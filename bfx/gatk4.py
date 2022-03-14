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
import logging
import re
import os

# MUGQIC Modules
import core
from core.config import global_config_parser
from core.job import Job
from . import gatk
from . import picard2

#################
#GATK4 - Read Data Manipulation

# only in GATK3
def realigner_target_creator(
    input,
    output,
    output_dir=[],
    input2=[],
    intervals=[],
    exclude_intervals=[]
    ):

    return gatk.realigner_target_creator(
        input,
        output,
        output_dir,
        input2,
        intervals,
        exclude_intervals
    )

# only in GATK3
def indel_realigner(
    input,
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

    return gatk.indel_realigner(
        input,
        target_intervals,
        input2,
        output,
        output_dir,
        output_norm_dep,
        output_tum_dep,
        intervals,
        exclude_intervals,
        optional
    )

def split_n_cigar_reads(
    input,
    output,
    intervals=[],
    exclude_intervals=[],
    interval_list=None
    ):

    if global_config_parser.param('gatk_split_N_trim', 'module_gatk').split("/")[2] < "4":
        return gatk.split_n_cigar_reads(
            input,
            output,
            intervals,
            exclude_intervals,
            interval_list
        )
    
    else:
        return Job(
            [input],
            [output],
            [
                ['gatk_split_N_trim', 'module_java'],
                ['gatk_split_N_trim', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  SplitNCigarReads {other_options} \\
  --TMP_DIR {tmp_dir} \\
  --reference_sequence {reference_sequence} \\
  --input_file {input} \\
  --out {output}{intervals}{exclude_intervals}""".format(
            tmp_dir=global_config_parser.param('gatk_split_N_trim', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_split_N_trim', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_split_N_trim', 'ram'),
            other_options=global_config_parser.param('gatk_split_N_trim', 'other_options', required=False),
            reference_sequence=global_config_parser.param('gatk_split_N_trim', 'reference', param_type='filepath'),
            input=input,
            output=output,
            intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
            interval_list=" \\\n --interval-padding 100 --intervals " + interval_list if interval_list else "",
            exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals)
        )
    )

def mark_duplicates(inputs, output, metrics_file, remove_duplicates="false"):
    if not isinstance(inputs, list):
        inputs = [inputs]

    return Job(
        inputs,
        [output, re.sub("\.([sb])am$", ".\\1ai", output), metrics_file],
        [
            ['gatk_mark_duplicates', 'module_java'],
            ['gatk_mark_duplicates', 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 MarkDuplicatesSpark \\
 --remove-all-duplicates {remove_duplicates} --read-validation-stringency SILENT --create-output-bam-index true \\
 --tmp-dir {tmp_dir} \\
 {inputs} \\
 --metrics-file {metrics_file} \\
 --spark-master local[{threads}] \\
 --output {output}""".format(
            tmp_dir=global_config_parser.param('gatk_mark_duplicates', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_mark_duplicates', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_mark_duplicates', 'ram'),
            remove_duplicates=remove_duplicates,
            inputs=" \\\n  ".join("--input " + input for input in inputs),
            threads=global_config_parser.param('gatk_mark_duplicates', 'threads', param_type='int'),
            output=output,
            metrics_file=metrics_file,
            max_records_in_ram=global_config_parser.param('gatk_mark_duplicates', 'max_records_in_ram', param_type='int')
        ),
        removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
    )

def base_recalibrator(input, output, intervals=None):
    if global_config_parser.param('gatk_base_recalibrator', 'module_gatk').split("/")[2] < "4":
        return gatk.base_recalibrator(input, output, intervals)
    else:
        return Job(
        [input, intervals],
        [output],
        [
            ['gatk_base_recalibrator', 'module_java'],
            ['gatk_base_recalibrator', 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  BaseRecalibratorSpark {options} \\
  --input {input} \\
  --reference {reference_sequence} {intervals} \\
  --known-sites {known_dbsnp} \\
  --known-sites {known_gnomad} \\
  --known-sites {known_mills} \\
  --spark-master local[{threads}] \\
  --output {output}""".format(
                tmp_dir=global_config_parser.param('gatk_base_recalibrator', 'tmp_dir'),
                java_other_options=global_config_parser.param('gatk_base_recalibrator', 'gatk4_java_options'),
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

def apply_bqsr(
    input,
    output,
    base_quality_score_recalibration
    ):

    if global_config_parser.param('gatk_apply_bqsr', 'module_gatk').split("/")[2] < "4":
        return gatk.print_reads(
            input,
            output,
            base_quality_score_recalibration
        )
    else:
        return Job(
            [
                input,
                base_quality_score_recalibration
            ],
            [output, re.sub(".bam", ".bam.bai", output)],
            [
                ['gatk_apply_bqsr', 'module_java'],
                ['gatk_apply_bqsr', 'module_gatk']
            ],
            command="""\
rm -rf {output}* && \\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  ApplyBQSRSpark {options} --create-output-bam-index true \\
  --input {input} \\
  --bqsr-recal-file {bqsr_file} \\
  --spark-master local[{threads}] \\
  --output {output}""".format(
            tmp_dir=global_config_parser.param('gatk_apply_bqsr', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_apply_bqsr', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_apply_bqsr', 'ram'),
            options=global_config_parser.param('gatk_apply_bqsr', 'options'),
            threads=global_config_parser.param('gatk_apply_bqsr', 'threads', param_type='int'),
            input=input,
            bqsr_file=base_quality_score_recalibration,
            output=output,
            )
        )

def print_reads(
    input,
    output,
    base_quality_score_recalibration
    ):

    if global_config_parser.param('print_reads', 'module_gatk').split("/")[2] < "4":
        return gatk.print_reads(
            input,
            output,
            base_quality_score_recalibration
        )
    else:
        return apply_bqsr(
            input,
            output,
            base_quality_score_recalibration
        )

#################
# GATK4 - Variant Manipulation

def cat_variants(
    variants,
    output=None
    ):

    if not isinstance(variants, list):
        variants = [variants]

    if global_config_parser.param('gatk_merge_vcfs', 'module_gatk').split("/")[2] < "4":
        return picard2.mergeVcfs(
            variants,
            output
        )
    else:
        return Job(
            variants,
            [output],
            [
                ['gatk_merge_vcfs', 'module_java'],
                ['gatk_merge_vcfs', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  MergeVcfs {options} \\
  --TMP_DIR {tmp_dir} \\
  --REFERENCE_SEQUENCE {reference}{variants} \\
  --OUTPUT {output}""".format(
                tmp_dir=global_config_parser.param('gatk_merge_vcfs', 'tmp_dir'),
                java_other_options=global_config_parser.param('gatk_merge_vcfs', 'gatk4_java_options'),
                ram=global_config_parser.param('gatk_merge_vcfs', 'ram'),
                options=global_config_parser.param('gatk_merge_vcfs', 'options'),
                reference=global_config_parser.param('gatk_merge_vcfs', 'genome_fasta', param_type='filepath'),
                variants="".join(" \\\n  --INPUT " + variant for variant in variants),
                output=output
            )
        )

def merge_stats(
        stats,
        output=None
        ):
        return Job(
            stats,
            [output],
            [
                ['gatk_merge_stats', 'module_java'],
                ['gatk_merge_stats', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  MergeMutectStats {options} \\
  {stats} \\
  --output {output}""".format(
                tmp_dir=global_config_parser.param('gatk_merge_stats', 'tmp_dir'),
                java_other_options=global_config_parser.param('gatk_merge_stats', 'gatk4_java_options'),
                ram=global_config_parser.param('gatk_merge_stats', 'ram'),
                options=global_config_parser.param('gatk_merge_stats', 'options'),
                stats="".join(" \\\n  --stats " + stat for stat in stats),
                output=output
            )
        )

def cnn_score_variants(
        input,
        output,
        input_bam,
        interval_list=None
        ):
    return Job(
        [input, input_bam],
        [output],
        [
            ['gatk_cnn_score_variants', 'module_java'],
            ['gatk_cnn_score_variants', 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CNNScoreVariants {options} \\
  --reference {reference_sequence} \\
  --variant {input} {input_bam} {interval_list} \\
  --output {output}""".format(
            tmp_dir=global_config_parser.param('gatk_cnn_score_variantss', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_cnn_score_variants', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_cnn_score_variants', 'ram'),
            options=global_config_parser.param('gatk_cnn_score_variants', 'options'),
            reference_sequence=global_config_parser.param('gatk_cnn_score_variants', 'reference', param_type='filepath'),
            input=input,
            input_bam=" \\\n  --input " + input_bam if input_bam else "",
            interval_list=" \\\n --intervals " + interval_list if interval_list else "",
            output=output
        )
    )

def filter_variant_tranches(
        input,
        output,
        interval_list=None
        ):
    return Job(
        [input],
        [output],
        [
            ['gatk_filter_variant_tranches', 'module_java'],
            ['gatk_filter_variant_tranches', 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  FilterVariantTranches {options} \\
  --resource {resource_hapmap} \\
  --resource {resource_mills} {interval_list} \\
  --variant {input} \\
  --output {output}""".format(
            tmp_dir=global_config_parser.param('gatk_filter_variant_tranches', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_filter_variant_tranches', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_filter_variant_tranches', 'ram'),
            options=global_config_parser.param('gatk_filter_variant_tranches', 'options'),
            resource_hapmap=global_config_parser.param('gatk_cnn_score_variants', 'hapmap', param_type='filepath'),
            resource_mills=global_config_parser.param('gatk_cnn_score_variants', 'mills', param_type='filepath'),
            input=input,
            interval_list=" \\\n --intervals " + interval_list if interval_list else "",
            output=output
        )
    )

# only in GATK3
def combine_variants(
    variants,
    output
    ):

    return gatk.combine_variants(
        variants,
        output
    )

def variant_annotator(
    input_normal,
    input_tumor,
    input_variants,
    output,
    intervals=[],
    exclude_intervals=[]
    ):

    return gatk.variant_annotator(
        input_normal,
        input_tumor,
        input_variants,
        output,
        intervals,
        exclude_intervals
    )

#####################
# GATK4 - Diagnostics and Quality Control

# only in GATK3
def depth_of_coverage(
    input,
    output_prefix,
    intervals
    ):

    return gatk.depth_of_coverage(
        input,
        output_prefix, 
        intervals
    )

def callable_loci(
    input,
    output,
    summary
    ):

    return gatk.callable_loci(
        input,
        output,
        summary
    )


#####################
# GATK4 - Short Variant Discovery

def haplotype_caller(
    inputs,
    output,
    intervals=[],
    exclude_intervals=[],
    interval_list=None
    ):

    interval_padding = config.param('gatk_haplotype_caller', 'interval_padding')

#added interval_padding as a varibale. Because in chipseq we don't need to add any padding to the peaks
    if not isinstance(inputs, list):
        inputs = [inputs]

    # Added this to check intervel_list (peak file) availability in the chip-seq pipeline
    inputs_list = inputs.copy()
    if not interval_list is None:
       inputs_list.extend([interval_list])

    if global_config_parser.param('gatk_haplotype_caller', 'module_gatk').split("/")[2] < "4":
        return gatk.haplotype_caller(
            inputs,
            output,
            intervals=intervals,
            exclude_intervals=exclude_intervals,
            interval_list=interval_list
        )
    else:
        return Job(
            #to track all files as input files replaced input with input_lists
            inputs_list,
            [output, output + ".tbi"],
            [
                ['gatk_haplotype_caller', 'module_java'],
                ['gatk_haplotype_caller', 'module_gatk']
            ],
            command="""\
gatk --java-options "{java_other_options} -Xmx{ram}" \\
  HaplotypeCaller {options} --native-pair-hmm-threads {threads} \\
  --reference {reference_sequence} \\
  --input {input} \\
  --output {output}{interval_padding} {interval_list}{intervals}{exclude_intervals}""".format(
                tmp_dir=global_config_parser.param('gatk_haplotype_caller', 'tmp_dir'),
                java_other_options=global_config_parser.param('gatk_haplotype_caller', 'gatk4_java_options'),
                ram=global_config_parser.param('gatk_haplotype_caller', 'ram'),
                options=global_config_parser.param('gatk_haplotype_caller', 'options'),
                threads=global_config_parser.param('gatk_haplotype_caller', 'threads'),
                reference_sequence=global_config_parser.param('gatk_haplotype_caller', 'genome_fasta', param_type='filepath'),
                interval_list=" --intervals " + interval_list if interval_list else "",
                interval_padding=" \\\n --interval-padding " + str(interval_padding)  if interval_padding else "",
                input=" \\\n  ".join(input for input in inputs),
                output=output,
                intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
                exclude_intervals="".join(
                    " \\\n  --exclude-intervals " + exclude_interval for exclude_interval in exclude_intervals)
            )
        )

def combine_gvcf(inputs, output, intervals=[], exclude_intervals=[]):
    if not isinstance(inputs, list):
        inputs = [inputs]
        
    if global_config_parser.param('gatk_haplotype_caller', 'module_gatk').split("/")[2] < "4":
        return gatk.combine_gvcf(inputs, output, intervals, exclude_intervals)
    else:
        
        return Job(
            inputs,
            [output],
            [
                ['gatk_combine_gvcf', 'module_java'],
                ['gatk_combine_gvcf', 'module_gatk']
            ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CombineGVCFs {other_options} \\
  --reference {reference_sequence} \\
  {input} \\
  --output {output}{intervals}{exclude_intervals}""".format(
                tmp_dir=global_config_parser.param('gatk_combine_gvcf', 'tmp_dir'),
                java_other_options=global_config_parser.param('gatk_combine_gvcf', 'gatk4_java_options'),
                ram=global_config_parser.param('gatk_combine_gvcf', 'ram'),
                other_options=global_config_parser.param('gatk_combine_gvcf', 'other_options', required=False),
                reference_sequence=global_config_parser.param('gatk_combine_gvcf', 'genome_fasta', param_type='filepath'),
                input="".join(" \\\n  --variant " + input for input in inputs),
                output=output,
                intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
                exclude_intervals="".join(
                    " \\\n  --exclude-intervals " + exclude_interval for exclude_interval in exclude_intervals)
        ))

def GenomicsDBImport(
    inputs,
    output,
    intervals=[],
    exclude_intervals=[]
    ):

    if not isinstance(inputs, list):
        inputs = [inputs]

    else:
        return Job(
            inputs,
            [output],
            [
                ['gatk_GenomicsDBImport', 'module_java'],
                ['gatk_GenomicsDBImport', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  GenomicsDBImport {other_options} \\
  {input} \\
  --output {output}{intervals}""".format(
                tmp_dir=global_config_parser.param('gatk_GenomicsDBImport', 'tmp_dir'),
                java_other_options=global_config_parser.param('gatk_GenomicsDBImport', 'gatk4_java_options'),
                ram=global_config_parser.param('gatk_GenomicsDBImport', 'ram'),
                other_options=global_config_parser.param('gatk_GenomicsDBImport', 'other_options', required=False),
                input="".join(" \\\n  --variant " + input for input in inputs),
                output=output,
                intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        ))

def genotype_gvcf(
    variants,
    output,
    options
    ):

    if not isinstance(variants, list):
        variants = [variants]

    if global_config_parser.param('gatk_genotype_gvcf', 'module_gatk').split("/")[2] < "4":
        return gatk.genotype_gvcf(
            variants,
            output,
            options
        )
    else:
        return Job(
            variants,
            [output],
            [
                ['gatk_genotype_gvcf', 'module_java'],
                ['gatk_genotype_gvcf', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  GenotypeGVCFs {options} \\
  --reference {reference_sequence}{variants} \\
  --output {output}""".format(
            tmp_dir=global_config_parser.param('gatk_genotype_gvcf', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_genotype_gvcf', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_genotype_gvcf', 'ram'),
            options=options,
            reference_sequence=global_config_parser.param('gatk_genotype_gvcf', 'genome_fasta', param_type='filepath'),
            variants="".join(" \\\n  --variant " + variant for variant in variants),
            output=output
        )
    )

def mutect2(inputNormal,
            normal_name,
            inputTumor,
            tumor_name,
            outputVCF,
            read_orientation,
            intervals=[],
            exclude_intervals=[],
            interval_list=None,
            ):

    if interval_list:
        inputs = [inputNormal, inputTumor, interval_list]

    else:
        inputs = [inputNormal, inputTumor]

    if global_config_parser.param('gatk_mutect2', 'module_gatk').split("/")[2] < "4":
        return gatk.mutect2(inputNormal,
                            normal_name,
                            inputTumor,
                            tumor_name,
                            outputVCF,
                            intervals,
                            exclude_intervals,
                            interval_list)
    else:
        return Job(
            inputs,
            [outputVCF, outputVCF + ".stats", read_orientation],
            [
                ['gatk_mutect2', 'module_java'],
                ['gatk_mutect2', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  Mutect2 {options} \\
  --f1r2-tar-gz {read_orientation} \\
  --reference {reference_sequence} \\
  --input {inputTumor} \\
  --tumor-sample {tumor_name} \\
  --input {inputNormal} \\
  --normal-sample {normal_name} \\
  --germline-resource {known_sites} \\
  --output {outputVCF}{interval_list}{intervals}{exclude_intervals}""".format(
        tmp_dir=global_config_parser.param('gatk_mutect', 'tmp_dir'),
        java_other_options=global_config_parser.param('gatk_mutect2', 'gatk4_java_options'),
        ram=global_config_parser.param('gatk_mutect2', 'ram'),
        options=global_config_parser.param('gatk_mutect2', 'options'),
        reference_sequence=global_config_parser.param('gatk_mutect2', 'genome_fasta', param_type='filepath'),
        read_orientation=read_orientation,
        known_sites=global_config_parser.param('gatk_mutect2', 'known_sites', param_type='filepath'),
        inputNormal=inputNormal,
        normal_name=normal_name,
        inputTumor=inputTumor,
        tumor_name=tumor_name,
        outputVCF=outputVCF,
        interval_list=" \\\n  --interval-padding 100 --intervals " + interval_list if interval_list else "",
        intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
        #pon=" --panel-of-normals " + global_config_parser.param('gatk_mutect2', 'pon', type='filepath') if global_config_parser.param('gatk_mutect2', 'pon', param_type='filepath') else "",
        exclude_intervals="".join(
            " \\\n  --exclude-intervals " + exclude_interval for exclude_interval in exclude_intervals)
        )
    )

#####################
# GATK4 - Variant Filtering

def get_pileup_summaries(
    input_bam,
    output,
    ):

    return Job(
        [input_bam],
        [output],
        [
            ['gatk_get_pileup_summaries', 'module_java'],
            ['gatk_get_pileup_summaries', 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  GetPileupSummaries {options} \\
  --input {input_bam} \\
  --variant {variants} --intervals {intervals} \\
  --output {output}""".format(
            tmp_dir=global_config_parser.param('gatk_get_pileup_summaries', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_get_pileup_summaries', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_get_pileup_summaries', 'ram'),
            options=global_config_parser.param('gatk_get_pileup_summaries', 'options'),
            variants=global_config_parser.param('gatk_get_pileup_summaries', 'known_sites', param_type='filepath'),
            intervals=global_config_parser.param('gatk_get_pileup_summaries', 'known_intervals', param_type='filepath'),
            input_bam=input_bam,
            output=output
        )
    )

def calculate_contamination(
        input,
        output,
        match_normal=None,
        tumor_segment=None
    ):

    return Job(
        [input, match_normal],
        [output, tumor_segment],
        [
            ['gatk_calculate_contamination', 'module_java'],
            ['gatk_calculate_contamination', 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CalculateContamination {options} \\
  --input {input} \\
  --output {output} \\
  {normal}{tumor_segment}""".format(
            tmp_dir=global_config_parser.param('gatk_calculate_contamination', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_calculate_contamination', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_calculate_contamination', 'ram'),
            options=global_config_parser.param('gatk_calculate_contamination', 'options'),
            input=input,
            normal=" \\\n --matched-normal " + match_normal if match_normal else "",
            tumor_segment=" \\\n --tumor-segmentation " + tumor_segment if tumor_segment else "",
            output=output,
        )
    )

def learn_read_orientation_model(
    inputs,
    output,
    ):

    return Job(
        inputs,
        [output],
        [
            ['gatk_learn_read_orientation_model', 'module_java'],
            ['gatk_learn_read_orientation_model', 'module_gatk']
        ],
    command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  LearnReadOrientationModel {options} \\
  {input} \\
  --output {output}""".format(
                tmp_dir=global_config_parser.param('gatk_learn_read_orientation_model', 'tmp_dir'),
                java_other_options=global_config_parser.param('gatk_learn_read_orientation_model', 'gatk4_java_options'),
                ram=global_config_parser.param('gatk_learn_read_orientation_model', 'ram'),
                options=global_config_parser.param('gatk_learn_read_orientation_model', 'options'),
                input="".join(" \\\n  --input " + input for input in inputs),
                output=output
            )
        )

def filter_mutect_calls(
    variants,
    output,
    contamination=None,
    segment=None,
    read_orientation=None
    ):

    return Job(
        [variants],
        [output],
        [
            ['gatk_filter_mutect_calls', 'module_java'],
            ['gatk_filter_mutect_calls', 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  FilterMutectCalls {options} \\
  --reference {reference} \\
  --variant {variants} \\
  {contamination_table} {segment_table} {read_orientation_model} \\
  --output {output}""".format(
                tmp_dir=global_config_parser.param('gatk_filter_mutect_calls', 'tmp_dir'),
                java_other_options=global_config_parser.param('gatk_filter_mutect_calls', 'gatk4_java_options'),
                ram=global_config_parser.param('gatk_filter_mutect_calls', 'ram'),
                options=global_config_parser.param('gatk_filter_mutect_calls', 'options'),
                reference=global_config_parser.param('gatk_filter_mutect_calls', 'genome_fasta', param_type='filepath'),
                variants=variants,
                contamination_table="  \\\n --contamination-table " + contamination if contamination else "",
                segment_table="  \\\n --tumor-segmentation " + segment if segment else "",
                read_orientation_model="  \\\n --ob-priors " + read_orientation if read_orientation else "",
                output=output
            )
        )

def variant_recalibrator(variants,
                         other_options,
                         recal_output,
                         tranches_output,
                         R_output,
                         small_sample_check=False):

    if not isinstance(variants, list):
        variants = [variants]

    if global_config_parser.param('gatk_variant_recalibrator', 'module_gatk').split("/")[2] < "4":
        return gatk.variant_recalibrator(variants, other_options, recal_output, tranches_output,
                                         R_output, small_sample_check=small_sample_check)
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
            [recal_output, tranches_output],
            [
                ['gatk_variant_recalibrator', 'module_java'],
                ['gatk_variant_recalibrator', 'module_gatk'],
                ['gatk_variant_recalibrator', 'module_R']
            ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  VariantRecalibrator {options} \\
  --reference {reference_sequence}{variants} \\
  {other_options} \\
  --output {recal_output} \\
  --tranches-file {tranches_output}""".format(
            tmp_dir=global_config_parser.param('gatk_variant_recalibrator', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_variant_recalibrator', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_variant_recalibrator', 'ram'),
            options=global_config_parser.param('gatk_variant_recalibrator', 'options'),
            reference_sequence=global_config_parser.param('gatk_variant_recalibrator', 'genome_fasta', param_type='filepath'),
            variants="".join(" \\\n  --variant " + variant for variant in variants),
            other_options=other_options,
            recal_output=recal_output,
            tranches_output=tranches_output,
            R_output=R_output
        ),
        removable_files=[
            recal_output,
            tranches_output,
            R_output
        ]
    )

def apply_recalibration(
    variants,
    recal_input,
    tranches_input,
    other_options,
    apply_recal_output
    ):

    if global_config_parser.param('gatk_apply_recalibration', 'module_gatk').split("/")[2] < "4":
        return gatk.apply_recalibration(
            variants,
            recal_input,
            tranches_input,
            other_options,
            apply_recal_output
        )
    else:
        return Job(
            [
                variants,
                recal_input,
                tranches_input
            ],
            [apply_recal_output],
            [
                ['gatk_apply_recalibration', 'module_java'],
                ['gatk_apply_recalibration', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  ApplyVQSR {options} \\
  --reference {reference_sequence} \\
  --variant {variants} \\
  {other_options} \\
  --tranches-file {tranches_input} \\
  --recal-file {recal_input} \\
  --output {output}""".format(
                tmp_dir=global_config_parser.param('gatk_apply_recalibration', 'tmp_dir'),
                java_other_options=global_config_parser.param('gatk_apply_recalibration', 'gatk4_java_options'),
                ram=global_config_parser.param('gatk_apply_recalibration', 'ram'),
                options=global_config_parser.param('gatk_apply_recalibration', 'options'),
                reference_sequence=global_config_parser.param('gatk_apply_recalibration', 'genome_fasta', param_type='filepath'),
                variants=variants,
                other_options=other_options,
                recal_input=recal_input,
                tranches_input=tranches_input,
                output=apply_recal_output
        ))

def variant_filtration(
    input,
    output,
    other_options
    ):

    return gatk.variant_filtration(
        input,
        output,
        other_options
    )


#####################
#  Copy Number Variant Discovery

def preprocessIntervals(input, output, intervals):
    return Job(
        [input],
        [output],
        [
            ['gatk_processIntervals', 'module_java'],
            ['gatk_processIntervals', 'module_gatk4']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  PreprocessIntervals {options} --interval-merging-rule OVERLAPPING_ONLY \\
  --reference {reference_sequence} {intervals} \\
  --bin-length {bin_length} \\
  --padding {padding} \\
  --output {output}""".format(
            tmp_dir=global_config_parser.param('gatk_processIntervals', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_processIntervals', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_processIntervals', 'ram'),
            options=global_config_parser.param('gatk_processIntervals', 'options'),
            reference_sequence=global_config_parser.param('gatk_processIntervals', 'genome_fasta', param_type='filepath'),
            intervals=" \\\n --intervals " + intervals if intervals else "",
            bin_length=global_config_parser.param('gatk_processIntervals', 'bin-length'),
            padding=global_config_parser.param('gatk_processIntervals', 'padding'),
            output=output
        )
    )

#####################
# PICARD imported functions

def build_bam_index(
    input,
    output
    ):

    if global_config_parser.param('build_bam_index', 'module_gatk').split("/")[2] < "4":
        return picard2.build_bam_index(
            input,
            output
        )
    else:
        return Job(
            [input],
            [output],
            [
                ['build_bam_index', 'module_java'],
                ['build_bam_index', 'module_gatk4']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  BuildBamIndex \\
  VALIDATION_STRINGENCY SILENT \\
  --TMP_DIR {tmp_dir} \\
  --INPUT {input} \\
  --OUTPUT {output} """.format(
                tmp_dir=global_config_parser.param('build_bam_index', 'tmp_dir'),
                java_other_options=global_config_parser.param('build_bam_index', 'gatk4_java_options'),
                ram=global_config_parser.param('build_bam_index', 'ram'),
                input=input,
                output=output,
        ))

def calculate_hs_metrics(
    input,
    output,
    intervals,
    reference_sequence=None
    ):

    baits_intervals = ""
    baits_intervals = global_config_parser.param('picard_calculate_hs_metrics', 'baits_intervals', required=False)
    
    if global_config_parser.param('picard_calculate_hs_metrics', 'module_gatk').split("/")[2] < "4":
        return picard2.calculate_hs_metrics(
            input,
            output,
            intervals,
            reference_sequence
        )
    else:

        return Job(
            [
                input,
                intervals
            ],
            [output],
            [
                ['picard_calculate_hs_metrics', 'module_java'],
                ['picard_calculate_hs_metrics', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 CollectHsMetrics \\
 --TMP_DIR {tmp_dir} \\
 --INPUT {input} \\
 --OUTPUT {output} \\
 --BAIT_INTERVALS {baits} \\
 --TARGET_INTERVALS {intervals} \\
 --REFERENCE_SEQUENCE {reference_sequence}""".format(
                tmp_dir=global_config_parser.param('picard_calculate_hs_metrics', 'tmp_dir'),
                java_other_options=global_config_parser.param('picard_calculate_hs_metrics', 'gatk4_java_options'),
                ram=global_config_parser.param('picard_calculate_hs_metrics', 'ram'),
                input=input,
                output=output,
                intervals=intervals,
                baits=baits_intervals if baits_intervals != "" else intervals,
                reference_sequence=reference_sequence if reference_sequence else global_config_parser.param('picard_calculate_hs_metrics', 'genome_fasta', param_type='filepath')
        ))

def collect_multiple_metrics(
    input,
    output,
    reference_sequence=None,
    library_type="PAIRED_END"
    ):

    if library_type == "PAIRED_END":
        outputs = [
            output + ".base_distribution_by_cycle_metrics",
            output + ".base_distribution_by_cycle.pdf",
            output + ".alignment_summary_metrics",
            output + ".insert_size_histogram.pdf",
            output + ".insert_size_metrics",
            output + ".quality_by_cycle_metrics",
            output + ".quality_by_cycle.pdf",
            output + ".quality_distribution_metrics",
            output + ".quality_distribution.pdf"
        ]
    else:
        outputs = [
            output + ".quality_by_cycle.pdf",
            output + ".alignment_summary_metrics",
            output + ".quality_by_cycle_metrics",
            output + ".quality_distribution_metrics",
            output + ".quality_distribution.pdf"
        ]

    if global_config_parser.param('picard_collect_multiple_metrics', 'module_gatk').split("/")[2] < "4":
        return picard2.collect_multiple_metrics(
            input,
            output,
            reference_sequence,
            library_type
        )
    
    else:
        
        return Job(
            [input],
            outputs,
            [
                ['picard_collect_multiple_metrics', 'module_java'],
                ['picard_collect_multiple_metrics', 'module_gatk'],
                ['picard_collect_multiple_metrics', 'module_R']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CollectMultipleMetrics \\
  --PROGRAM CollectAlignmentSummaryMetrics \\
  --PROGRAM CollectInsertSizeMetrics \\
  --VALIDATION_STRINGENCY SILENT \\
  --TMP_DIR {tmp_dir} \\
  --REFERENCE_SEQUENCE {reference_sequence} \\
  --INPUT {input} \\
  --OUTPUT {output} \\
  --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                tmp_dir=global_config_parser.param('picard_collect_multiple_metrics', 'tmp_dir'),
                java_other_options=global_config_parser.param('picard_collect_multiple_metrics', 'gatk4_java_options'),
                ram=global_config_parser.param('picard_collect_multiple_metrics', 'ram'),
                reference_sequence=reference_sequence if reference_sequence else global_config_parser.param('picard_collect_multiple_metrics', 'genome_fasta', param_type='filepath'),
                input=input,
                output=output,
                max_records_in_ram=global_config_parser.param('picard_collect_multiple_metrics', 'max_records_in_ram', param_type='int')
        ))

def collect_sequencing_artifacts_metrics(
    input,
    output,
    annotation_flat=None,
    reference_sequence=None
    ):

    output_dep = output + ".bait_bias_summary_metrics"
    
    if global_config_parser.param('picard_collect_sequencing_artifacts_metrics', 'module_gatk').split("/")[2] < "4":
        return picard2.collect_sequencing_artifacts_metrics(
            input,
            output,
            annotation_flat,
            reference_sequence
        )
    
    else:
        return Job(
            [input],
            [output_dep],
            [
                ['picard_collect_sequencing_artifacts_metrics', 'module_java'],
                ['picard_collect_sequencing_artifacts_metrics', 'module_gatk'],
                ['picard_collect_sequencing_artifacts_metrics', 'module_R']
            ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 CollectSequencingArtifactMetrics \\
  --VALIDATION_STRINGENCY SILENT {options} \\
  --INPUT {input} \\
  --OUTPUT {output} \\
  --REFERENCE_SEQUENCE {reference} \\
  --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
            options=global_config_parser.param('picard_collect_sequencing_artifacts_metrics', 'options'),
            tmp_dir=global_config_parser.param('picard_collect_sequencing_artifacts_metrics', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_collect_sequencing_artifacts_metrics', 'gatk4_java_options'),
            ram=global_config_parser.param('picard_collect_sequencing_artifacts_metrics', 'ram'),
            input=input,
            output=output,
            reference=reference_sequence if reference_sequence else global_config_parser.param('picard_collect_sequencing_artifacts_metrics', 'genome_fasta'),
            max_records_in_ram=global_config_parser.param('picard_collect_sequencing_artifacts_metrics', 'max_records_in_ram', param_type='int')
    ))

def convert_sequencing_artifacts_metrics(
    input,
    output,
    annotation_flat=None,
    reference_sequence=None
    ):

    input_dep = input + ".bait_bias_summary_metrics"
    
    if global_config_parser.param('picard_convert_sequencing_artifacts_metrics', 'module_gatk').split("/")[2] < "4":
        return picard2.convert_sequencing_artifacts_metrics(
            input,
            output,
            annotation_flat,
            reference_sequence
        )
    else:
        return Job(
            [input_dep],
            [output],
            [
                ['picard_convert_sequencing_artifacts_metrics', 'module_java'],
                ['picard_convert_sequencing_artifacts_metrics', 'module_gatk'],
                ['picard_convert_sequencing_artifacts_metrics', 'module_R']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  ConvertSequencingArtifactToOxoG \\
  --VALIDATION_STRINGENCY SILENT  \\
  --TMP_DIR {tmp_dir} \\
  --INPUT_BASE {input} \\
  --OUTPUT_BASE {output} \\
  --REFERENCE_SEQUENCE {reference}""".format(
                tmp_dir=global_config_parser.param('picard_convert_sequencing_artifacts_metrics', 'tmp_dir'),
                java_other_options=global_config_parser.param('picard_convert_sequencing_artifacts_metrics', 'gatk4_java_options'),
                ram=global_config_parser.param('picard_convert_sequencing_artifacts_metrics', 'ram'),
                input=input,
                output=output,
                reference=reference_sequence if reference_sequence else global_config_parser.param('picard_convert_sequencing_artifacts_metrics', 'genome_fasta'),
            )
        )
    
def collect_oxog_metrics(
    input,
    output,
    annotation_flat=None,
    reference_sequence=None
    ):
    
    if global_config_parser.param('picard_collect_oxog_metrics', 'module_gatk').split("/")[2] < "4":
        return picard2.collect_oxog_metrics(
            input,
            output,
            annotation_flat,
            reference_sequence
        )
    else:
        return Job(
            [input],
            [output],
            [
                ['picard_collect_oxog_metrics', 'module_java'],
                ['picard_collect_oxog_metrics', 'module_gatk'],
                ['picard_collect_oxog_metrics', 'module_R']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CollectOxoGMetrics \\
  --VALIDATION_STRINGENCY SILENT  \\
  --TMP_DIR {tmp_dir} \\
  --INPUT {input} \\
  --OUTPUT {output} \\
  {dbsnp} \\
  --REFERENCE_SEQUENCE {reference} \\
  --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                tmp_dir=global_config_parser.param('picard_collect_oxog_metrics', 'tmp_dir'),
                java_other_options=global_config_parser.param('picard_collect_oxog_metrics', 'gatk4_java_options'),
                ram=global_config_parser.param('picard_collect_oxog_metrics', 'ram'),
                input=input,
                output=output,
                dbsnp="--DB_SNP " + global_config_parser.param('picard_collect_oxog_metrics', 'known_variants') if global_config_parser.param('picard_collect_oxog_metrics', 'known_variants') else "",
                reference=reference_sequence if reference_sequence else global_config_parser.param('picard_collect_oxog_metrics', 'genome_fasta'),
                max_records_in_ram=global_config_parser.param('picard_collect_oxog_metrics', 'max_records_in_ram', param_type='int')
        ))

def collect_gcbias_metrics(
    input,
    output,
    chart,
    summary_file,
    annotation_flat=None,
    reference_sequence=None
    ):
    
    if global_config_parser.param('picard_collect_gcbias_metrics', 'module_gatk').split("/")[2] < "4":
        return picard2.collect_gcbias_metrics(
            input,
            output,
            chart,
            summary_file,
            annotation_flat,
            reference_sequence
        )
    else:
        return Job(
            [input],
            [output],
            [
                ['picard_collect_gcbias_metrics', 'module_java'],
                ['picard_collect_gcbias_metrics', 'module_gatk'],
                ['picard_collect_gcbias_metrics', 'module_R']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CollectGcBiasMetrics \\
  --VALIDATION_STRINGENCY SILENT \\
  --ALSO_IGNORE_DUPLICATES TRUE \\
  --TMP_DIR {tmp_dir} \\
  --INPUT {input} \\
  --OUTPUT {output} \\
  --CHART {chart} \\
  --SUMMARY_OUTPUT {summary_file} \\
  --REFERENCE_SEQUENCE {reference} \\
  --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
            tmp_dir=global_config_parser.param('picard_collect_gcbias_metrics', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_collect_gcbias_metrics', 'gatk4_java_options'),
            ram=global_config_parser.param('picard_collect_gcbias_metrics', 'ram'),
            input=input,
            output=output,
            chart=chart,
            summary_file=summary_file,
            reference=reference_sequence if reference_sequence else global_config_parser.param('picard_collect_gcbias_metrics',
                                                                                 'genome_fasta'),
            max_records_in_ram=global_config_parser.param('picard_collect_gcbias_metrics', 'max_records_in_ram', param_type='int')
        )
    )

def fix_mate_information(
    input,
    output
    ):
    
    if global_config_parser.param('picard_fix_mate_information', 'module_gatk').split("/")[2] < "4":
        return picard2.fix_mate_information(
            input,
            output
        )
    else:
        return Job(
            [input],
            [output],
            [
                ['picard_fix_mate_information', 'module_java'],
                ['picard_fix_mate_information', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 FixMateInformation \\
 --VALIDATION_STRINGENCY SILENT CREATE_INDEX true SORT_ORDER coordinate \\
 --TMP_DIR {tmp_dir} \\
 --INPUT {input} \\
 --OUTPUT {output} \\
 --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                tmp_dir=global_config_parser.param('picard_fix_mate_information', 'tmp_dir'),
                java_other_options=global_config_parser.param('picard_fix_mate_information', 'gatk4_java_options'),
                ram=global_config_parser.param('picard_fix_mate_information', 'ram'),
                input=input,
                output=output,
                max_records_in_ram=global_config_parser.param('picard_fix_mate_information', 'max_records_in_ram', param_type='int')
            ),
            removable_files=[
                output,
                re.sub("\.([sb])am$", ".\\1ai", output),
                output + ".md5"
            ]
        )

def picard_mark_duplicates(
    inputs,
    output,
    metrics_file,
    remove_duplicates="false"
    ):

    if not isinstance(inputs, list):
        inputs = [inputs]
        
    if global_config_parser.param('picard_mark_duplicates', 'module_gatk').split("/")[2] < "4":
        return picard2.mark_duplicates(
            inputs,
            output,
            metrics_file,
            remove_duplicates
        )
    else:
        return Job(
            inputs,
            [
                output,
                re.sub("\.([sb])am$", ".\\1ai", output),
                metrics_file
            ],
            [
                ['picard_mark_duplicates', 'module_java'],
                ['picard_mark_duplicates', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 MarkDuplicates \\
 --REMOVE_DUPLICATES {remove_duplicates} \\
 --VALIDATION_STRINGENCY SILENT \\
 --CREATE_INDEX true \\
 --TMP_DIR {tmp_dir} \\
 {inputs} \\
 --OUTPUT {output} \\
 --METRICS_FILE {metrics_file} \\
 --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                tmp_dir=global_config_parser.param('picard_mark_duplicates', 'tmp_dir'),
                java_other_options=global_config_parser.param('picard_mark_duplicates', 'gatk4_java_options'),
                ram=global_config_parser.param('picard_mark_duplicates', 'ram'),
                remove_duplicates=remove_duplicates,
                inputs=" \\\n  ".join("--INPUT " + input for input in inputs),
                output=output,
                metrics_file=metrics_file,
                max_records_in_ram=global_config_parser.param('picard_mark_duplicates', 'max_records_in_ram', param_type='int')
            ),
            removable_files=[
                output,
                re.sub("\.([sb])am$", ".\\1ai", output),
                output + ".md5"
            ]
        )

def mark_duplicates_mate_cigar(inputs,
                           output,
                           metrics_file,
                           remove_duplicates="false"):
    if not isinstance(inputs, list):
        inputs = [inputs]
    
    if global_config_parser.param('mark_duplicates_mate_cigar', 'module_gatk').split("/")[2] < "4":
        return picard2.mark_duplicates_mate_cigar(inputs,
                                       output,
                                       metrics_file,
                                       remove_duplicates)
    else:
        return Job(
            inputs,
            [output, re.sub("\.([sb])am$", ".\\1ai", output), metrics_file],
            [
                ['mark_duplicates_mate_cigar', 'module_java'],
                ['mark_duplicates_mate_cigar', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 MarkDuplicatesWithMateCigar \\
 --REMOVE_DUPLICATES {remove_duplicates} \\
 --VALIDATION_STRINGENCY SILENT \\
 --CREATE_INDEX true \\
 --TMP_DIR {tmp_dir} \\
 {inputs} \\
 --OUTPUT {output} \\
 --METRICS_FILE {metrics_file} \\
 --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                tmp_dir=global_config_parser.param('mark_duplicates_mate_cigar', 'tmp_dir'),
                java_other_options=global_config_parser.param('mark_duplicates_mate_cigar', 'java_other_options'),
                ram=global_config_parser.param('mark_duplicates_mate_cigar', 'ram'),
                remove_duplicates=remove_duplicates,
                inputs=" \\\n  ".join("--INPUT " + input for input in inputs),
                output=output,
                metrics_file=metrics_file,
                max_records_in_ram=global_config_parser.param('mark_duplicates_mate_cigar', 'max_records_in_ram', param_type='int')
            ),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
        )

def picard_mark_duplicates_mate_cigar(inputs,
                           output,
                           metrics_file,
                           remove_duplicates="false"):
    if not isinstance(inputs, list):
        inputs = [inputs]
    
    if global_config_parser.param('picard_mark_duplicates_mate_cigar', 'module_gatk').split("/")[2] < "4":
        return picard2.mark_duplicates_mate_cigar(inputs,
                                       output,
                                       metrics_file,
                                       remove_duplicates)
    else:
        return Job(
            inputs,
            [output, re.sub("\.([sb])am$", ".\\1ai", output), metrics_file],
            [
                ['picard_mark_duplicates_mate_cigar', 'module_java'],
                ['picard_mark_duplicates_mate_cigar', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 MarkDuplicatesWithMateCigar \\
 --REMOVE_DUPLICATES {remove_duplicates} \\
 --VALIDATION_STRINGENCY SILENT \\
 --CREATE_INDEX true \\
 --TMP_DIR {tmp_dir} \\
 {inputs} \\
 --OUTPUT {output} \\
 --METRICS_FILE {metrics_file} \\
 --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                tmp_dir=global_config_parser.param('picard_mark_duplicates_mate_cigar', 'tmp_dir'),
                java_other_options=global_config_parser.param('picard_mark_duplicates_mate_cigar', 'gatk4_java_options'),
                ram=global_config_parser.param('picard_mark_duplicates_mate_cigar', 'ram'),
                remove_duplicates=remove_duplicates,
                inputs=" \\\n  ".join("--INPUT " + input for input in inputs),
                output=output,
                metrics_file=metrics_file,
                max_records_in_ram=global_config_parser.param('picard_mark_duplicates_mate_cigar', 'max_records_in_ram', param_type='int')
            ),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
        )

def merge_sam_files(inputs,
                    output):
    
    if not isinstance(inputs, list):
        inputs = [inputs]

    if global_config_parser.param('picard_merge_sam_files', 'module_gatk').split("/")[2] < "4":
        return picard2.merge_sam_files(
            inputs,
            output
        )
    else:
        
        return Job(
            inputs,
            [
                output,
                re.sub("\.([sb])am$", ".\\1ai", output)
            ],
            [
                ['picard_merge_sam_files', 'module_java'],
                ['picard_merge_sam_files', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 MergeSamFiles \\
 --VALIDATION_STRINGENCY SILENT \\
 --ASSUME_SORTED true \\
 --CREATE_INDEX true \\
 --TMP_DIR {tmp_dir} \\
 {inputs} \\
 --OUTPUT {output} \\
 --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                tmp_dir=global_config_parser.param('picard_merge_sam_files', 'tmp_dir'),
                java_other_options=global_config_parser.param('picard_merge_sam_files', 'gatk4_java_options'),
                ram=global_config_parser.param('picard_merge_sam_files', 'ram'),
                inputs=" \\\n ".join(["--INPUT " + input for input in inputs]),
                output=output,
                max_records_in_ram=global_config_parser.param('picard_merge_sam_files', 'max_records_in_ram', param_type='int')
            ),
            removable_files=[
                output,
                re.sub("\.([sb])am$", ".\\1ai", output)
            ]
        )

# Reorder BAM/SAM files based on reference/dictionary
def reorder_sam(
    input,
    output
    ):

    return Job(
        [input],
        [output],
        [
            ['picard_reorder_sam', 'module_java'],
            ['picard_reorder_sam', 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 ReorderSam \\
 --VALIDATION_STRINGENCY SILENT \\
 --CREATE_INDEX true \\
 --TMP_DIR {tmp_dir} \\
 --INPUT {input} \\
 --OUTPUT {output} \\
 --REFERENCE {reference} \\
 --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                tmp_dir=global_config_parser.param('picard_reorder_sam', 'tmp_dir'),
                java_other_options=global_config_parser.param('picard_reorder_sam', 'gatk4_java_options'),
                ram=global_config_parser.param('picard_reorder_sam', 'ram'),
                input=input,
                output=output,
                reference=global_config_parser.param('picard_reorder_sam', 'genome_fasta', param_type='filepath'),
                max_records_in_ram=global_config_parser.param('picard_reorder_sam', 'max_records_in_ram', param_type='int')
            ),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output)]
        )

# Convert SAM/BAM file to fastq format
def sam_to_fastq(
    input,
    fastq,
    second_end_fastq=None
    ):

    return Job(
        [input],
        [
            fastq,
            second_end_fastq
        ],
        [
            ['picard_sam_to_fastq', 'module_java'],
            ['picard_sam_to_fastq', 'module_gatk']
        ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 SamToFastq \\
 --VALIDATION_STRINGENCY LENIENT \\
 --CREATE_MD5_FILE TRUE \\
 --INPUT {input} \\
 --FASTQ {fastq}{second_end_fastq}""".format(
                tmp_dir=global_config_parser.param('picard_sam_to_fastq', 'tmp_dir'),
                java_other_options=global_config_parser.param('picard_sam_to_fastq', 'gatk4_java_options'),
                ram=global_config_parser.param('picard_sam_to_fastq', 'ram'),
                input=input,
                fastq=fastq,
                second_end_fastq=" \\\n --SECOND_END_FASTQ " + second_end_fastq if second_end_fastq else ""
            ),
            removable_files=[
                fastq,
                second_end_fastq
            ]
        )

def sort_sam(
    input,
    output,
    sort_order="coordinate",
    ini_section='picard_sort_sam'
    ):
    
    if global_config_parser.param(ini_section, 'module_gatk').split("/")[2] < "4":
        return picard2.sort_sam(
            input,
            output,
            sort_order,
            ini_section
        )
    else:
        return Job(
            [input],
            # Add SAM/BAM index as output only when writing a coordinate-sorted BAM file
            [
                output,
                re.sub("\.([sb])am$", ".\\1ai", output) if sort_order == "coordinate" else None
            ],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  SortSam \\
  --VALIDATION_STRINGENCY SILENT \\
  --CREATE_INDEX true \\
  --TMP_DIR {tmp_dir} \\
  --INPUT {input} \\
  --OUTPUT {output} \\
  --SORT_ORDER {sort_order} \\
  --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                tmp_dir=global_config_parser.param(ini_section, 'tmp_dir'),
                java_other_options=global_config_parser.param(ini_section, 'gatk4_java_options'),
                ram=global_config_parser.param(ini_section, 'ram'),
                input=input,
                output=output,
                sort_order=sort_order,
                max_records_in_ram=global_config_parser.param(ini_section, 'max_records_in_ram', param_type='int')
            ),
            removable_files=[
                output,
                re.sub("\.([sb])am$", ".\\1ai", output) if sort_order == "coordinate" else None
            ]
        )

def sort_vcfs(
    inputs,
    output,
    ini_section='picard_sort_vcf'
    ):

    if not isinstance(inputs, list):
        inputs = [inputs]
        
    return Job(
        inputs,
        # Add SAM/BAM index as output only when writing a coordinate-sorted BAM file
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 SortVcf \\
 VALIDATION_STRINGENCY SILENT \\
 TMP_DIR {tmp_dir} \\
 {inputs} \\
 OUTPUT {output} \\
 SEQUENCE_DICTIONARY {seq_dict}""".format(
                tmp_dir=global_config_parser.param(ini_section, 'tmp_dir'),
                java_other_options=global_config_parser.param(ini_section, 'gatk4_java_options'),
                ram=global_config_parser.param(ini_section, 'ram'),
                inputs=" \\\n  ".join(["INPUT " + input for input in inputs]),
                output=output,
                seq_dict=global_config_parser.param(ini_section, 'genome_dictionary', param_type='filepath')
            )
        )

def collect_rna_metrics(
    input,
    output,
    annotation_flat=None,
    reference_sequence=None
    ):

    return Job(
        [input],
        # collect specific RNA metrics (exon rate, strand specificity, etc...)
        [output],
        [
            ['picard_collect_rna_metrics', 'module_java'],
            ['picard_collect_rna_metrics', 'module_gatk'],
            ['picard_collect_rna_metrics', 'module_R']
        ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 CollectRnaSeqMetrics \\
 --VALIDATION_STRINGENCY SILENT  \\
 --TMP_DIR {tmp_dir} \\
 --INPUT {input} \\
 --OUTPUT {output} \\
 --REF_FLAT {ref_flat} \\
 --STRAND_SPECIFICITY {strand_specificity} \\
 --MINIMUM_LENGTH {min_length} \\
 --REFERENCE_SEQUENCE {reference} \\
 --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                tmp_dir=global_config_parser.param('picard_collect_rna_metrics', 'tmp_dir'),
                java_other_options=global_config_parser.param('picard_collect_rna_metrics', 'gatk4_java_options'),
                ram=global_config_parser.param('picard_collect_rna_metrics', 'ram'),
                input=input,
                output=output,
                ref_flat=annotation_flat if annotation_flat else global_config_parser.param('picard_collect_rna_metrics', 'annotation_flat'),
                strand_specificity=global_config_parser.param('picard_collect_rna_metrics', 'strand_info'),
                min_length=global_config_parser.param('picard_collect_rna_metrics', 'minimum_length', param_type='int'),
                reference=reference_sequence if reference_sequence else global_config_parser.param('picard_collect_rna_metrics', 'genome_fasta'),
                max_records_in_ram=global_config_parser.param('picard_collect_rna_metrics', 'max_records_in_ram', param_type='int')
            )
        )


def crosscheck_fingerprint(inputs,
                           output):

    matrix=output + ".matrix"
        
    return Job(
        inputs,
        [output],
        [
            ['gatk_crosscheck_fingerprint', 'module_java'],
            ['gatk_crosscheck_fingerprint', 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CrosscheckFingerprints {options} \\
  --VALIDATION_STRINGENCY SILENT \\
  --TMP_DIR {tmp_dir} \\
  {inputs} \\
  --HAPLOTYPE_MAP {haplotype_database} \\
  --LOD_THRESHOLD {lod_threshold} \\
  --OUTPUT {output} \\
  --MATRIX_OUTPUT {matrix}""".format(
                tmp_dir=global_config_parser.param('gatk_crosscheck_fingerprint', 'tmp_dir'),
                options=global_config_parser.param('gatk_crosscheck_fingerprint', 'options'),
                java_other_options=global_config_parser.param('gatk_crosscheck_fingerprint', 'gatk4_java_options'),
                haplotype_database=global_config_parser.param('gatk_crosscheck_fingerprint', 'haplotype_database'),
                lod_threshold=global_config_parser.param('gatk_crosscheck_fingerprint', 'lod_threshold'),
                ram=global_config_parser.param('gatk_crosscheck_fingerprint', 'ram'),
                inputs=" \\\n  ".join("--INPUT " + str(input) for input in inputs),
                output=output,
                matrix=matrix,
            )
        )

def cluster_crosscheck_metrics(input,
                               output):
        return Job(
            [input],
            [output],
            [
                ['gatk_cluster_crosscheck_metrics', 'module_java'],
                ['gatk_cluster_crosscheck_metrics', 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  ClusterCrosscheckMetrics {options} \\
  --VALIDATION_STRINGENCY SILENT \\
  --TMP_DIR {tmp_dir} \\
  --INPUT {input} \\
  --LOD_THRESHOLD {lod_threshold} \\
  --OUTPUT {output} """.format(
                tmp_dir=global_config_parser.param('gatk_cluster_crosscheck_metrics', 'tmp_dir'),
                options=global_config_parser.param('gatk_cluster_crosscheck_metrics', 'options'),
                java_other_options=global_config_parser.param('gatk_cluster_crosscheck_metrics', 'gatk4_java_options'),
                lod_threshold=global_config_parser.param('gatk_cluster_crosscheck_metrics', 'lod_threshold'),
                ram=global_config_parser.param('gatk_cluster_crosscheck_metrics', 'ram'),
                input=input,
                output=output,
            )
        )

def bed2interval_list(
    dictionary,
    bed,
    output
    ):
	
    if global_config_parser.param('gatk_bed2interval_list', 'module_gatk').split("/")[2] < "4":
        return gatk.bed2interval_list(
            dictionary,
            bed,
            output
            )

    return Job(
        [dictionary, bed],
        [output],
        [
            ['gatk_bed2interval_list', 'module_java'],
            ['gatk_bed2interval_list', 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  BedToIntervalList \\
  --INPUT {bed} \\
  --SEQUENCE_DICTIONARY {dictionary} \\
  --OUTPUT {output}""".format(
            tmp_dir=global_config_parser.param('gatk_bed2interval_list', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_bed2interval_list', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_bed2interval_list', 'ram'),
            dictionary=dictionary if dictionary else global_config_parser.param('gatk_bed2interval_list', 'genome_dictionary', param_type='filepath'),
            bed=bed,
            output=output
            )
        )

def interval_list2bed(input, output):
    return Job(
        [input],
        [output],
        [
            ['gatk_interval_list2bed', 'module_java'],
            ['gatk_interval_list2bed', 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  IntervalListToBed \\
  --INPUT {input} \\
  --OUTPUT {output}""".format(
            tmp_dir=global_config_parser.param('gatk_interval_list2bed', 'tmp_dir'),
            java_other_options=global_config_parser.param('gatk_interval_list2bed', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_interval_list2bed', 'ram'),
            input=input,
            output=output
            )
        )

def preProcessInterval(reference,
                       intervals,
                       output,
                       options = None):
#                  exclude_intervals=None):

    return Job(
        [intervals],
        [output],
        [
            ['gatk_preProcessInterval', 'module_java'],
            ['gatk_preProcessInterval', 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  PreprocessIntervals {options} \\
  --reference {reference} \\
  --intervals {intervals} \\
  --output {output}""".format(
            tmp_dir=global_config_parser.param('gatk_preProcessInterval', 'tmp_dir'),
            options=options if options else global_config_parser.param('gatk_preProcessInterval', 'options'),
            java_other_options=global_config_parser.param('gatk_preProcessInterval', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_preProcessInterval', 'ram'),
            reference=reference if reference else global_config_parser.param('gatk_preProcessInterval', 'genome_fasta', param_type='filepath'),
            intervals=intervals,
#            exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals),
            output=output
        )
    )

def splitInterval(intervals,
                  output,
                  jobs,
                  options = None):
#                  exclude_intervals=None):

    interval_list = []
    for idx in range(jobs):
     interval_list.append(
         os.path.join(output, str(idx).zfill(4) + "-scattered.interval_list")
         )

    return Job(
        [intervals],
        interval_list,
        [
            ['gatk_splitInterval', 'module_java'],
            ['gatk_splitInterval', 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  SplitIntervals {options} \\
  --scatter-count {jobs} \\
  --reference {reference} \\
  --intervals {intervals} \\
  --output {output}""".format(
            tmp_dir=global_config_parser.param('gatk_splitInterval', 'tmp_dir'),
            options=options if options else global_config_parser.param('gatk_splitInterval', 'options'),
            jobs=jobs,
            java_other_options=global_config_parser.param('gatk_splitInterval', 'gatk4_java_options'),
            ram=global_config_parser.param('gatk_splitInterval', 'ram'),
            reference=global_config_parser.param('gatk_splitInterval', 'genome_fasta', param_type='filepath'),
            intervals=intervals,
#            exclude_intervals="".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals),
            output=output
        )
    )
