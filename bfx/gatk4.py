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
import gatk
import picard2

#################
#GATK4 - Read Data Manipulation

# only in GATK3
def realigner_target_creator(input, output, input2=[], intervals=[], exclude_intervals=[]):
	return gatk.realigner_target_creator(input, output, input2, intervals, exclude_intervals)

# only in GATK3
def indel_realigner(input, target_intervals, input2=[], output=[], output_norm_dep=[], output_tum_dep=[], intervals=[],
                    exclude_intervals=[], optional=[]):
	return gatk.indel_realigner(input, target_intervals, input2, output, output_norm_dep, output_tum_dep, intervals,
                    exclude_intervals, optional)

def split_n_cigar_reads(input, output, intervals=[], exclude_intervals=[], interval_list=None):
	if config.param('gatk_split_N_trim', 'module_gatk').split("/")[2] < "4":
		return gatk.split_n_cigar_reads(input, output, intervals, exclude_intervals, interval_list)
	
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

def mark_duplicates(inputs, output, metrics_file, remove_duplicates="false"):

	if not isinstance(inputs, list):
		inputs=[inputs]
 
	return Job(
        inputs,
        [output, re.sub("\.([sb])am$", ".\\1ai", output), metrics_file],
        [
            ['gatk_mark_duplicates', 'module_java'],
            ['gatk_mark_duplicates', 'module_gatk']
        ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  MarkDuplicatesGATK \\
  --REMOVE_DUPLICATES={remove_duplicates} --VALIDATION_STRINGENCY SILENT --CREATE_INDEX true \\
  --TMP_DIR {tmp_dir} \\
  {inputs} \\
  --metrics-file {metrics_file} \\
  --output {output} \\
  --MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
            tmp_dir=config.param('gatk_mark_duplicates', 'tmp_dir'),
            java_other_options=config.param('gatk_mark_duplicates', 'java_other_options'),
            ram=config.param('gatk_mark_duplicates', 'ram'),
            remove_duplicates=remove_duplicates,
            inputs=" \\\n  ".join("--input " + input for input in inputs),
            output=output,
            metrics_file=metrics_file,
            max_records_in_ram=config.param('gatk_mark_duplicates', 'max_records_in_ram', type='int')
            ),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
        )

def base_recalibrator(input, output, intervals=None):
	
	if config.param('base_recalibrator', 'module_gatk').split("/")[2] < "4":
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

def apply_bqsr(input, output, base_quality_score_recalibration):
	if config.param('gatk_apply_bqsr', 'module_gatk').split("/")[2] < "4":
		return gatk.print_reads(input, output, base_quality_score_recalibration)
	else:
		return Job(
			[input, base_quality_score_recalibration],
			[output],
			[
				['gatk_apply_bqsr', 'module_java'],
				['gatk_apply_bqsr', 'module_gatk']
			],
		command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  ApplyBQSRSpark {options} \\
  --reference {reference_sequence} \\
  --input {input} \\
  --bqsr-recal-file {bqsr_file} \\
  --spark-master local[{threads}] \\
  --output {output}""".format(
			tmp_dir=config.param('gatk_apply_bqsr', 'tmp_dir'),
			java_other_options=config.param('gatk_apply_bqsr', 'java_other_options'),
			ram=config.param('gatk_apply_bqsr', 'ram'),
			options=config.param('gatk_apply_bqsr', 'options'),
			threads=config.param('gatk_apply_bqsr', 'threads', type='int'),
			reference_sequence=config.param('gatk_apply_bqsr', 'genome_2bit', type='filepath'),
			input=input,
			bqsr_file=base_quality_score_recalibration,
			output=output,
		)
	)

def print_reads(input, output, base_quality_score_recalibration):
	if config.param('print_reads', 'module_gatk').split("/")[2] < "4":
		return gatk.print_reads(input, output, base_quality_score_recalibration)
	else:
		return apply_bqsr(input, output, base_quality_score_recalibration)

#################
# GATK4 - Variant Manipulation

def cat_variants(variants, output=None):
	if config.param('gatk_cat_variants', 'module_gatk').split("/")[2] < "4":
		return gatk.cat_variants(variants, output)
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
				tmp_dir=config.param('gatk_merge_vcfs', 'tmp_dir'),
				java_other_options=config.param('gatk_merge_vcfs', 'java_other_options'),
				ram=config.param('gatk_merge_vcfs', 'ram'),
				options=config.param('gatk_merge_vcfs', 'options'),
				reference=config.param('gatk_merge_vcfs', 'genome_fasta', type='filepath'),
				variants="".join(" \\\n  --INPUT " + variant for variant in variants),
				output=output
			)
		)

# only in GATK3
def combine_variants(variants, output):
	return gatk.combine_variants(variants, output)

def variant_annotator(input_normal, input_tumor, input_variants, output, intervals=[], exclude_intervals=[]):
	return gatk.variant_annotator(input_normal, input_tumor, input_variants, output, intervals, exclude_intervals)


#####################
# GATK4 - Diagnostics and Quality Control

# only in GATK3
def depth_of_coverage(input, output_prefix, intervals):
	return gatk.depth_of_coverage(input, output_prefix, intervals)

def callable_loci(input, output, summary):
	return gatk.callable_loci(input, output, summary)


#####################
# GATK4 - Short Variant Discovery

def haplotype_caller(inputs, output, intervals=[], exclude_intervals=[], interval_list=None):
	if not isinstance(inputs, list):
		inputs = [inputs, interval_list]
	
	if config.param('gatk_haplotype_caller', 'module_gatk').split("/")[2] < "4":
		return gatk.haplotype_caller(inputs, output, intervals, exclude_intervals, interval_list)
	else:
		return Job(
			inputs,
			[output],
			[
				['gatk_haplotype_caller', 'module_java'],
				['gatk_haplotype_caller', 'module_gatk']
			],
			command="""\
gatk --java-options "{java_other_options} -Xmx{ram}" \\
  HaplotypeCaller {options} \\
  --reference {reference_sequence} \\
  --input {input} \\
  --output {output}{intervals}{exclude_intervals}""".format(
				tmp_dir=config.param('gatk_haplotype_caller', 'tmp_dir'),
				java_other_options=config.param('gatk_haplotype_caller', 'java_other_options'),
				ram=config.param('gatk_haplotype_caller', 'ram'),
				options=config.param('gatk_haplotype_caller', 'options'),
				reference_sequence=config.param('gatk_haplotype_caller', 'genome_fasta', type='filepath'),
				input=" \\\n  ".join(input for input in inputs),
				output=output,
				interval_list=" \\\n --interval-padding 100 --intervals " + interval_list if interval_list else "",
				intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
				exclude_intervals="".join(
					" \\\n  --exclude-intervals " + exclude_interval for exclude_interval in exclude_intervals)
			)
		)

def combine_gvcf(inputs, output, intervals=[], exclude_intervals=[]):
	if not isinstance(inputs, list):
		inputs = [inputs]
	if config.param('gatk_haplotype_caller', 'module_gatk').split("/")[2] < "4":
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
			tmp_dir=config.param('gatk_combine_gvcf', 'tmp_dir'),
			java_other_options=config.param('gatk_combine_gvcf', 'java_other_options'),
			ram=config.param('gatk_combine_gvcf', 'ram'),
			other_options=config.param('gatk_combine_gvcf', 'other_options', required=False),
			reference_sequence=config.param('gatk_combine_gvcf', 'genome_fasta', type='filepath'),
			input="".join(" \\\n  --variant " + input for input in inputs),
			output=output,
			intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
			exclude_intervals="".join(
				" \\\n  --exclude-intervals " + exclude_interval for exclude_interval in exclude_intervals)
		)
	)

def genotype_gvcf(variants, output, options):
	if config.param('genotype_gvcf', 'module_gatk').split("/")[2] < "4":
		return gatk.genotype_gvcf(variants, output, options)
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
			tmp_dir=config.param('gatk_genotype_gvcf', 'tmp_dir'),
			java_other_options=config.param('gatk_genotype_gvcf', 'java_other_options'),
			ram=config.param('gatk_genotype_gvcf', 'ram'),
			options=options,
			reference_sequence=config.param('gatk_genotype_gvcf', 'genome_fasta', type='filepath'),
			variants="".join(" \\\n  --variant " + variant for variant in variants),
			output=output
		)
	)


def mutect2(inputNormal, normal_name, inputTumor, tumor_name, outputVCF, intervals=[], exclude_intervals=[], interval_list=None):
	if config.param('gatk_mutect2', 'module_gatk').split("/")[2] < "4":
		return gatk.mutect2(inputNormal, normal_name, inputTumor, tumor_name, outputVCF, intervals, exclude_intervals, interval_list)
	else:
		return Job(
			[inputNormal, inputTumor, interval_list],
			[outputVCF],
			[
				['gatk_mutect2', 'module_java'],
				['gatk_mutect2', 'module_gatk']
			],
		command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  Mutect2 {options} \\
  --reference {reference_sequence} \\
  --input {inputTumor} \\
  --tumor-sample {tumor_name} \\
  --input {inputNormal} \\
  --normal-sample {normal_name} \\
  --germline-resource {known_sites} \\
  --output {outputVCF}{interval_list}{intervals}{exclude_intervals}""".format(
			tmp_dir=config.param('gatk_mutect', 'tmp_dir'),
			java_other_options=config.param('gatk_mutect2', 'java_other_options'),
			ram=config.param('gatk_mutect2', 'ram'),
			options=config.param('gatk_mutect2', 'options'),
			reference_sequence=config.param('gatk_mutect2', 'genome_fasta', type='filepath'),
			known_sites=config.param('gatk_mutect2', 'known_sites', type='filepath'),
			inputNormal=inputNormal,
			normal_name=normal_name,
			inputTumor=inputTumor,
			tumor_name=tumor_name,
			outputVCF=outputVCF,
			interval_list=" \\\n  --interval-padding 100 --intervals " + interval_list if interval_list else "",
			intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
			exclude_intervals="".join(
				" \\\n  --exclude-intervals " + exclude_interval for exclude_interval in exclude_intervals)
		)
	)


#####################
# GATK4 - Variant Filtering

def get_pileup_summaries(input_bam, output):
	return Job(
		[input_bam],
		[output],
		[
			['get_pileup_summaries', 'module_java'],
			['get_pileup_summaries', 'module_gatk']
		],
		command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  GetPileupSummaries {options} \\
  --input {input_bam} \\
  --variant {variants} \\
  --intervals {intervals} \\
  --output {output}""".format(
			tmp_dir=config.param('get_pileup_summaries', 'tmp_dir'),
			java_other_options=config.param('get_pileup_summaries', 'java_other_options'),
			ram=config.param('get_pileup_summaries', 'ram'),
			options=config.param('get_pileup_summaries', 'options'),
			variants=config.param('get_pileup_summaries', 'known_sites', type='filepath'),
			intervals=config.param('get_pileup_summaries', 'intervals', type='filepath'),
			input_bam=input_bam,
			output=output,
		)
	)

def calculate_contamination(input, output, match_normal=None):
	return Job(
		[input],
		[output],
		[
			['get_pileup_summaries', 'module_java'],
			['get_pileup_summaries', 'module_gatk']
		],
		command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CalculateContamination {options} \\
  --tmp_dir {tmp_dir} \\
  --input {input} \\
  {normal} \\
  --output {output}""".format(
			tmp_dir=config.param('get_pileup_summaries', 'tmp_dir'),
			java_other_options=config.param('get_pileup_summaries', 'java_other_options'),
			ram=config.param('get_pileup_summaries', 'ram'),
			options=config.param('get_pileup_summaries', 'options'),
			input=input,
			normal=" \\\n --matched-normal " + match_normal if match_normal else "",
			output=output,
		)
	)

def filter_mutect_calls(variants, output, contamination=None):
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
  --tmp_dir {tmp_dir} \\
  --variant {variants} \\
  {contamination_table} \\
  --output {output}""".format(
				tmp_dir=config.param('gatk_filter_mutect_call', 'tmp_dir'),
				java_other_options=config.param('gatk_filter_mutect_call', 'java_other_options'),
				ram=config.param('gatk_filter_mutect_call', 'ram'),
				options=config.param('gatk_filter_mutect_call', 'options'),
				variants=variants,
				contamination_table=" \\\n --contamination-table " + contamination if contamination else "",
				output=output
			)
		)

def variant_recalibrator(variants, other_options, recal_output, tranches_output, R_output):
	if config.param('gatk_variant_recalibrator', 'module_gatk').split("/")[2] < "4":
		return gatk.variant_recalibrator(variants, other_options, recal_output, tranches_output, R_output)
	else:
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
			tmp_dir=config.param('gatk_variant_recalibrator', 'tmp_dir'),
			java_other_options=config.param('gatk_variant_recalibrator', 'java_other_options'),
			ram=config.param('gatk_variant_recalibrator', 'ram'),
			options=config.param('gatk_variant_recalibrator', 'options'),
			reference_sequence=config.param('gatk_variant_recalibrator', 'genome_fasta', type='filepath'),
			variants="".join(" \\\n  --variant " + variant for variant in variants),
			other_options=other_options,
			#tmp_dir="--TMP_DIR " + tmp_dir if tmp_dir else "",
			recal_output=recal_output,
			tranches_output=tranches_output,
			R_output=R_output
		),
		removable_files=[recal_output, tranches_output, R_output]
	)


def apply_recalibration(variants, recal_input, tranches_input, other_options, apply_recal_output):
	if config.param('gatk_variant_recalibrator', 'module_gatk').split("/")[2] < "4":
		return gatk.apply_recalibration(variants, recal_input, tranches_input, other_options, apply_recal_output)
	else:
		return Job(
			[variants, recal_input, tranches_input],
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

def variant_filtration(input, output, other_options):
	return gatk.variant_filtration(input, output, other_options)


#####################
# PICARD imported functions

def build_bam_index(input, output):
	if config.param('build_bam_index', 'module_gatk').split("/")[2] < "4":
		return picard2.build_bam_index(input, output)
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
  VALIDATION_STRINGENCY=SILENT \\
  --TMP_DIR {tmp_dir} \\
  --INPUT={input} \\
  --OUTPUT={output} """.format(
				tmp_dir=config.param('build_bam_index', 'tmp_dir'),
				java_other_options=config.param('build_bam_index', 'java_other_options'),
				ram=config.param('build_bam_index', 'ram'),
				input=input,
				output=output,
			)
		)


def calculate_hs_metrics(input, output, intervals, reference_sequence=None):
	baits_intervals = ""
	baits_intervals = config.param('picard_calculate_hs_metrics', 'baits_intervals', required=False)
	
	if config.param('picard_calculate_hs_metrics', 'module_gatk').split("/")[2] < "4":
		return picard2.calculate_hs_metrics(input, output, intervals, reference_sequence)
	else:

		return Job(
			[input, intervals],
			[output],
			[
				['picard_calculate_hs_metrics', 'module_java'],
				['picard_calculate_hs_metrics', 'module_gatk']
			],
			command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CollectHsMetrics \\
  --TMP_DIR={tmp_dir} \\
  --INPUT={input} \\
  --OUTPUT={output} \\
  --BAIT_INTERVALS={baits} \\
  --TARGET_INTERVALS={intervals} \\
  --REFERENCE_SEQUENCE={reference_sequence}""".format(
				tmp_dir=config.param('picard_calculate_hs_metrics', 'tmp_dir'),
				java_other_options=config.param('picard_calculate_hs_metrics', 'java_other_options'),
				ram=config.param('picard_calculate_hs_metrics', 'ram'),
				input=input,
				output=output,
				intervals=intervals,
				baits=baits_intervals if baits_intervals != "" else intervals,
				reference_sequence=reference_sequence if reference_sequence else config.param(
					'picard_calculate_hs_metrics', 'genome_fasta', type='filepath')
			)
		)


def collect_multiple_metrics(input, output, reference_sequence=None, library_type="PAIRED_END"):
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

	if config.param('picard_collect_multiple_metrics', 'module_gatk').split("/")[2] < "4":
		return picard2.collect_multiple_metrics(input, output, reference_sequence, library_type)
	
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
  --PROGRAM=CollectAlignmentSummaryMetrics \\
  --PROGRAM=CollectInsertSizeMetrics \\
  --VALIDATION_STRINGENCY=SILENT \\
  --TMP_DIR={tmp_dir} \\
  --REFERENCE_SEQUENCE={reference_sequence} \\
  --INPUT={input} \\
  --OUTPUT={output} \\
  --MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
				tmp_dir=config.param('picard_collect_multiple_metrics', 'tmp_dir'),
				java_other_options=config.param('picard_collect_multiple_metrics', 'java_other_options'),
				ram=config.param('picard_collect_multiple_metrics', 'ram'),
				reference_sequence=reference_sequence if reference_sequence else config.param(
					'picard_collect_multiple_metrics', 'genome_fasta', type='filepath'),
				input=input,
				output=output,
				max_records_in_ram=config.param('picard_collect_multiple_metrics', 'max_records_in_ram', type='int')
			)
		)


def collect_sequencing_artifacts_metrics(input, output, annotation_flat=None, reference_sequence=None):
	output_dep = output + ".bait_bias_summary_metrics"
	
	if config.param('picard_collect_sequencing_artifacts_metrics', 'module_gatk').split("/")[2] < "4":
		return picard2.collect_sequencing_artifacts_metrics(input, output, annotation_flat, reference_sequence)
	
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
  --VALIDATION_STRINGENCY=SILENT {options} \\
  --TMP_DIR={tmp_dir} \\
  --INPUT={input} \\
  --OUTPUT={output} \\
  --REFERENCE_SEQUENCE={reference} \\
  --MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
			options=config.param('picard_collect_sequencing_artifacts_metrics', 'options'),
			tmp_dir=config.param('picard_collect_sequencing_artifacts_metrics', 'tmp_dir'),
			java_other_options=config.param('picard_collect_sequencing_artifacts_metrics', 'java_other_options'),
			ram=config.param('picard_collect_sequencing_artifacts_metrics', 'ram'),
			input=input,
			output=output,
			reference=reference_sequence if reference_sequence else config.param(
				'picard_collect_sequencing_artifacts_metrics', 'genome_fasta'),
			max_records_in_ram=config.param('picard_collect_sequencing_artifacts_metrics', 'max_records_in_ram',
			                                type='int')
		)
	)


def convert_sequencing_artifacts_metrics(input, output, annotation_flat=None, reference_sequence=None):
	input_dep = input + ".bait_bias_summary_metrics"
	
	if config.param('picard_convert_sequencing_artifacts_metrics', 'module_gatk').split("/")[2] < "4":
		return picard2.convert_sequencing_artifacts_metrics(input, output, annotation_flat, reference_sequence)
	
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
  --VALIDATION_STRINGENCY=SILENT  \\
  --TMP_DIR={tmp_dir} \\
  --INPUT_BASE={input} \\
  --OUTPUT_BASE={output} \\
  --REFERENCE_SEQUENCE={reference}""".format(
				tmp_dir=config.param('picard_convert_sequencing_artifacts_metrics', 'tmp_dir'),
				java_other_options=config.param('picard_convert_sequencing_artifacts_metrics', 'java_other_options'),
				ram=config.param('picard_convert_sequencing_artifacts_metrics', 'ram'),
				input=input,
				output=output,
				reference=reference_sequence if reference_sequence else config.param(
					'picard_convert_sequencing_artifacts_metrics', 'genome_fasta'),
			)
		)


def collect_oxog_metrics(input, output, annotation_flat=None, reference_sequence=None):
	
	if config.param('picard_collect_oxog_metrics', 'module_gatk').split("/")[2] < "4":
		return picard2.collect_oxog_metrics(input, output, annotation_flat, reference_sequence)
	
	else:
		
		return Job(
			[input],
			[output],
			[
				['picard_collect_sequencing_artifacts_metrics', 'module_java'],
				['picard_collect_sequencing_artifacts_metrics', 'module_gatk'],
				['picard_collect_sequencing_artifacts_metrics', 'module_R']
			],
		command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CollectOxoGMetrics \\
  --VALIDATION_STRINGENCY=SILENT  \\
  --TMP_DIR={tmp_dir} \\
  --INPUT={input} \\
  --OUTPUT={output} \\
  --DB_SNP={dbsnp} \\
  --REFERENCE_SEQUENCE={reference} \\
  --MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
			tmp_dir=config.param('picard_collect_oxog_metrics', 'tmp_dir'),
			java_other_options=config.param('picard_collect_oxog_metrics', 'java_other_options'),
			ram=config.param('picard_collect_oxog_metrics', 'ram'),
			input=input,
			output=output,
			dbsnp=config.param('picard_collect_oxog_metrics', 'known_variants'),
			reference=reference_sequence if reference_sequence else config.param('picard_collect_oxog_metrics',
			                                                                     'genome_fasta'),
			max_records_in_ram=config.param('picard_collect_oxog_metrics', 'max_records_in_ram', type='int')
		)
	)


def collect_gcbias_metrics(input, output, chart, summary_file, annotation_flat=None, reference_sequence=None):
	
	if config.param('picard_collect_gcbias_metrics', 'module_gatk').split("/")[2] < "4":
		return picard2.collect_gcbias_metrics(input, output, chart, summary_file, annotation_flat, reference_sequence)
	
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
  --VALIDATION_STRINGENCY=SILENT \\
  --ALSO_IGNORE_DUPLICATES=TRUE \\
  --TMP_DIR={tmp_dir} \\
  --INPUT={input} \\
  --OUTPUT={output} \\
  --CHART={chart} \\
  --SUMMARY_OUTPUT={summary_file} \\
  --REFERENCE_SEQUENCE={reference} \\
  --MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
			tmp_dir=config.param('picard_collect_gcbias_metrics', 'tmp_dir'),
			java_other_options=config.param('picard_collect_gcbias_metrics', 'java_other_options'),
			ram=config.param('picard_collect_gcbias_metrics', 'ram'),
			input=input,
			output=output,
			chart=chart,
			summary_file=summary_file,
			reference=reference_sequence if reference_sequence else config.param('picard_collect_gcbias_metrics',
			                                                                     'genome_fasta'),
			max_records_in_ram=config.param('picard_collect_gcbias_metrics', 'max_records_in_ram', type='int')
		)
	)

def fix_mate_information(input, output):
	
	if config.param('fix_mate_information', 'module_gatk').split("/")[2] < "4":
		return picard2.fix_mate_information(input, output)
	else:
	
		return Job(
			[input],
			[output],
			[
				['fixmate', 'module_java'],
				['fixmate', 'module_gatk']
			],
			command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  FixMateInformation \\
  --VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate \\
  --TMP_DIR={tmp_dir} \\
  --INPUT={input} \\
  --OUTPUT={output} \\
  --MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
				tmp_dir=config.param('picard_fix_mate_information', 'tmp_dir'),
				java_other_options=config.param('picard_fix_mate_information', 'java_other_options'),
				ram=config.param('picard_fix_mate_information', 'ram'),
				input=input,
				output=output,
				max_records_in_ram=config.param('picard_fix_mate_information', 'max_records_in_ram', type='int')
			),
			removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
		)


def picard_mark_duplicates(inputs, output, metrics_file, remove_duplicates="false"):
	if not isinstance(inputs, list):
		inputs = [inputs]
		
	if config.param('picard_mark_duplicates', 'module_gatk').split("/")[2] < "4":
		return picard2.mark_duplicates(inputs, output, metrics_file, remove_duplicates)
	
	else:

		return Job(
		inputs,
		[output, re.sub("\.([sb])am$", ".\\1ai", output), metrics_file],
		[
			['picard_mark_duplicates', 'module_java'],
			['picard_mark_duplicates', 'module_gatk']
		],
			command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  MarkDuplicates \\
  --REMOVE_DUPLICATES={remove_duplicates} \\
  --VALIDATION_STRINGENCY=SILENT \\
  --CREATE_INDEX=true \\
  --TMP_DIR={tmp_dir} \\
  {inputs} \\
  --OUTPUT={output} \\
  --METRICS_FILE={metrics_file} \\
  --MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
				tmp_dir=config.param('picard_mark_duplicates', 'tmp_dir'),
				java_other_options=config.param('picard_mark_duplicates', 'java_other_options'),
				ram=config.param('picard_mark_duplicates', 'ram'),
				remove_duplicates=remove_duplicates,
				inputs=" \\\n  ".join("--INPUT=" + input for input in inputs),
				output=output,
				metrics_file=metrics_file,
				max_records_in_ram=config.param('picard_mark_duplicates', 'max_records_in_ram', type='int')
			),
			removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
		)


def merge_sam_files(inputs, output):
	if not isinstance(inputs, list):
		inputs = [inputs]

	if config.param('fix_mate_information', 'module_gatk').split("/")[2] < "4":
		return picard2.merge_sam_files(inputs, output)
	else:
		
		return Job(
			inputs,
			[output, re.sub("\.([sb])am$", ".\\1ai", output)],
			[
				['picard_merge_sam_files', 'module_java'],
				['picard_merge_sam_files', 'module_gatk']
			],
			command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  MergeSamFiles \\
  --VALIDATION_STRINGENCY=SILENT \\
  --ASSUME_SORTED=true \\
  --CREATE_INDEX=true \\
  --TMP_DIR={tmp_dir} \\
  {inputs} \\
  --OUTPUT={output} \\
  --MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
				tmp_dir=config.param('picard_merge_sam_files', 'tmp_dir'),
				java_other_options=config.param('picard_merge_sam_files', 'java_other_options'),
				ram=config.param('picard_merge_sam_files', 'ram'),
				inputs=" \\\n ".join(["--INPUT=" + input for input in inputs]),
				output=output,
				max_records_in_ram=config.param('picard_merge_sam_files', 'max_records_in_ram', type='int')
			),
			removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output)]
		)


# Reorder BAM/SAM files based on reference/dictionary
def reorder_sam(input, output):
	return Job(
		[input],
		[output],
		[
			['reorder_sam', 'module_java'],
			['reorder_sam', 'module_gatk']
		],
			command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  ReorderSam \\
  --VALIDATION_STRINGENCY=SILENT \\
  --CREATE_INDEX=true \\
  --TMP_DIR={tmp_dir} \\
  --INPUT={input} \\
  --OUTPUT={output} \\
  --REFERENCE={reference} \\
  --MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
				tmp_dir=config.param('picard_reorder_sam', 'tmp_dir'),
				java_other_options=config.param('picard_reorder_sam', 'java_other_options'),
				ram=config.param('picard_reorder_sam', 'ram'),
				input=input,
				output=output,
				reference=config.param('picard_reorder_sam', 'genome_fasta', type='filepath'),
				max_records_in_ram=config.param('picard_reorder_sam', 'max_records_in_ram', type='int')
			),
			removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output)]
		)


# Convert SAM/BAM file to fastq format
def sam_to_fastq(input, fastq, second_end_fastq=None):
	return Job(
		[input],
		[fastq, second_end_fastq],
		[
			['picard_sam_to_fastq', 'module_java'],
			['picard_sam_to_fastq', 'module_gatk']
		],
			command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  SamToFastq \\
  --VALIDATION_STRINGENCY=LENIENT \\
  --CREATE_MD5_FILE=TRUE \\
  --TMP_DIR={tmp_dir} \\
  --INPUT={input} \\
  --FASTQ={fastq}{second_end_fastq}""".format(
				tmp_dir=config.param('picard_sam_to_fastq', 'tmp_dir'),
				java_other_options=config.param('picard_sam_to_fastq', 'java_other_options'),
				ram=config.param('picard_sam_to_fastq', 'ram'),
				input=input,
				fastq=fastq,
				second_end_fastq=" \\\n --SECOND_END_FASTQ=" + second_end_fastq if second_end_fastq else ""
			),
			removable_files=[fastq, second_end_fastq]
		)


def sort_sam(input, output, sort_order="coordinate", ini_section='picard_sort_sam'):
	
	if config.param('sort_sam', 'module_gatk').split("/")[2] < "4":
		return picard2.sort_sam(input, output)
	
	else:

		return Job(
			[input],
			# Add SAM/BAM index as output only when writing a coordinate-sorted BAM file
			[output, re.sub("\.([sb])am$", ".\\1ai", output) if sort_order == "coordinate" else None],
			[
				[ini_section, 'module_java'],
				[ini_section, 'module_gatk']
			],
			command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  SortSam \\
  --VALIDATION_STRINGENCY=SILENT \\
  --CREATE_INDEX=true \\
  --TMP_DIR={tmp_dir} \\
  --INPUT={input} \\
  --OUTPUT={output} \\
  --SORT_ORDER={sort_order} \\
  --MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
				tmp_dir=config.param(ini_section, 'tmp_dir'),
				java_other_options=config.param(ini_section, 'java_other_options'),
				ram=config.param(ini_section, 'ram'),
				input=input,
				output=output,
				sort_order=sort_order,
				max_records_in_ram=config.param(ini_section, 'max_records_in_ram', type='int')
			),
			removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output) if sort_order == "coordinate" else None]
		)


def sort_vcfs(inputs, output, ini_section='picard_sort_vcf'):
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
  VALIDATION_STRINGENCY=SILENT \\
  TMP_DIR={tmp_dir} \\
  {inputs} \\
  OUTPUT={output} \\
  SEQUENCE_DICTIONARY={seq_dict}""".format(
				tmp_dir=config.param(ini_section, 'tmp_dir'),
				java_other_options=config.param(ini_section, 'java_other_options'),
				ram=config.param(ini_section, 'ram'),
				inputs=" \\\n  ".join(["INPUT=" + input for input in inputs]),
				output=output,
				seq_dict=config.param(ini_section, 'genome_dictionary', type='filepath')
			)
		)


def collect_rna_metrics(input, output, annotation_flat=None, reference_sequence=None):
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
  --VALIDATION_STRINGENCY=SILENT  \\
  --TMP_DIR={tmp_dir} \\
  --INPUT={input} \\
  --OUTPUT={output} \\
  --REF_FLAT={ref_flat} \\
  --STRAND_SPECIFICITY={strand_specificity} \\
  --MINIMUM_LENGTH={min_length} \\
  --REFERENCE_SEQUENCE={reference} \\
  --MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
				tmp_dir=config.param('picard_collect_rna_metrics', 'tmp_dir'),
				java_other_options=config.param('picard_collect_rna_metrics', 'java_other_options'),
				ram=config.param('picard_collect_rna_metrics', 'ram'),
				input=input,
				output=output,
				ref_flat=annotation_flat if annotation_flat else config.param('picard_collect_rna_metrics',
				                                                              'annotation_flat'),
				strand_specificity=config.param('picard_collect_rna_metrics', 'strand_info'),
				min_length=config.param('picard_collect_rna_metrics', 'minimum_length', type='int'),
				reference=reference_sequence if reference_sequence else config.param('picard_collect_rna_metrics',
				                                                                     'genome_fasta'),
				max_records_in_ram=config.param('picard_collect_rna_metrics', 'max_records_in_ram', type='int')
			)
		)

