################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes. If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import re
import os

# GenPipes Modules
from ..core.config import global_conf
from ..core.job import Job
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
    exclude_intervals=[],
    fix_encoding=[]
    ):

    return gatk.realigner_target_creator(
        input,
        output,
        output_dir,
        input2,
        intervals,
        exclude_intervals,
        fix_encoding
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
    optional=[],
    fix_encoding=[]
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
        optional,
        fix_encoding
    )

def split_n_cigar_reads(
    input,
    output,
    intervals=[],
    exclude_intervals=[],
    interval_list=None,
    ini_section='gatk_split_N_trim'
    ):

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
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
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  SplitNCigarReads {other_options} \\
  --tmp-dir {tmp_dir} \\
  --reference {reference_sequence} \\
  --input {input} \\
  --output {output}{intervals}{exclude_intervals}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                other_options=global_conf.global_get(ini_section, 'other_options', required=False),
                reference_sequence=global_conf.global_get(ini_section, 'reference', param_type='filepath'),
                input=input,
                output=output,
                intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
                interval_list=" \\\n --interval-padding 100 --intervals " + interval_list if interval_list else "",
                exclude_intervals="".join(" \\\n  --exclude-intervals " + exclude_interval for exclude_interval in exclude_intervals)
        )
    )

def mark_duplicates_spark(
    inputs,
    output,
    metrics_file,
    remove_duplicates="false",
    ini_section='gatk_mark_duplicates_spark'
    ):
    if not isinstance(inputs, list):
        inputs = [inputs]

    return Job(
        inputs,
        [output, re.sub(r"\.([sb])am$", ".\\1ai", output), metrics_file],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
rm -rf {output}.part && \\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 MarkDuplicatesSpark \\
 --remove-all-duplicates {remove_duplicates} --read-validation-stringency SILENT --create-output-bam-index true \\
 --tmp-dir {tmp_dir} \\
 {inputs} \\
 --metrics-file {metrics_file} \\
 --spark-master local[{threads}] \\
 --output {output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            remove_duplicates=remove_duplicates,
            inputs=" \\\n  ".join("--input " + input for input in inputs),
            threads=global_conf.global_get(ini_section, 'threads', param_type='int'),
            output=output,
            metrics_file=metrics_file,
            max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
        ),
    )

def base_recalibrator(
    input,
    output,
    intervals=None,
    ini_section='gatk_base_recalibrator'
    ):
    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return gatk.base_recalibrator(input, output, intervals)
    else:
        return Job(
            [input, intervals],
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
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
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                options=global_conf.global_get(ini_section, 'options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                threads=global_conf.global_get(ini_section, 'threads', param_type='int'),
                input=input,
                intervals=" \\\n  --intervals " + intervals if intervals else "",
                reference_sequence=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
                known_dbsnp=global_conf.global_get(ini_section, 'known_dbsnp', param_type='filepath'),
                known_gnomad=global_conf.global_get(ini_section, 'known_gnomad', param_type='filepath'),
                known_mills=global_conf.global_get(ini_section, 'known_mills', param_type='filepath'),
                output=output
            ),
        )

def apply_bqsr(
    input,
    output,
    base_quality_score_recalibration,
    ini_section='gatk_apply_bqsr'
    ):

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
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
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
            ],
            command="""\
rm -rf {output}* && \\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  ApplyBQSRSpark {options} --create-output-bam-index true \\
  --input {input} \\
  --bqsr-recal-file {bqsr_file} \\
  --spark-master local[{threads}] \\
  --output {output} && \\
chmod 660 {outputs}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                options=global_conf.global_get(ini_section, 'options'),
                threads=global_conf.global_get(ini_section, 'threads', param_type='int'),
                input=input,
                bqsr_file=base_quality_score_recalibration,
                output=output,
                outputs=f"{output} {re.sub('.bam', '.bam.bai', output)}"
            )
        )

def print_reads(
    input,
    output,
    base_quality_score_recalibration
    ):

    if global_conf.global_get('print_reads', 'module_gatk').split("/")[2] < "4":
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
    output=None,
    ini_section='gatk_merge_vcfs'
    ):

    if not isinstance(variants, list):
        variants = [variants]

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return gatk.cat_variants(
            variants,
            output
        )
    else:
        return Job(
            variants,
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  MergeVcfs {options} \\
  --TMP_DIR {tmp_dir} \\
  --REFERENCE_SEQUENCE {reference}{variants} \\
  --OUTPUT {output}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                options=global_conf.global_get(ini_section, 'options'),
                reference=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
                variants="".join(" \\\n  --INPUT " + variant for variant in variants),
                output=output
            )
        )

def merge_stats(
        stats,
        output=None,
        ini_section='gatk_merge_stats'
        ):
    return Job(
        stats,
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  MergeMutectStats {options} \\
  {stats} \\
  --output {output}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                options=global_conf.global_get(ini_section, 'options'),
                stats="".join(" \\\n  --stats " + stat for stat in stats),
                output=output
            )
        )

def cnn_score_variants(
        input,
        output,
        input_bam,
        interval_list=None,
        ini_section='gatk_cnn_score_variants'
        ):
    return Job(
        [input, input_bam],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CNNScoreVariants {options} \\
  --reference {reference_sequence} \\
  --variant {input} {input_bam} {interval_list} \\
  --output {output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            options=global_conf.global_get(ini_section, 'options'),
            reference_sequence=global_conf.global_get(ini_section, 'reference', param_type='filepath'),
            input=input,
            input_bam=" \\\n  --tensor-type read-tensor --input " + input_bam if input_bam else "",
            interval_list=" \\\n --intervals " + interval_list if interval_list else "",
            output=output
        )
    )

def filter_variant_tranches(
        input,
        output,
        interval_list=None,
        ini_section='gatk_filter_variant_tranches'
        ):
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  FilterVariantTranches {options} \\
  --resource {resource_hapmap} \\
  --resource {resource_mills} {interval_list} \\
  --variant {input} \\
  --output {output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            options=global_conf.global_get(ini_section, 'options'),
            resource_hapmap=global_conf.global_get(ini_section, 'hapmap', param_type='filepath'),
            resource_mills=global_conf.global_get(ini_section, 'mills', param_type='filepath'),
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
    input,
    output,
    interval_list,
    intervals=None,
    exclude_intervals=None,
    ini_section='gatk_haplotype_caller'
    ):

    interval_padding = global_conf.global_get(ini_section, 'interval_padding')

#added interval_padding as a variable. Because in chipseq we don't need to add any padding to the peaks
    inputs = []
    if not isinstance(input, list):
        inputs = [input]

    # Added this to check intervel_list (peak file) availability in the chip-seq pipeline
    inputs_list = inputs.copy()
    if not interval_list is None:
       inputs_list.append(interval_list)

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return gatk.haplotype_caller(
            input,
            output,
            interval_list
        )
    else:
        return Job(
            #to track all files as input files replaced input with input_lists
            inputs_list,
            [output, output + ".tbi"],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
            ],
            command="""\
gatk --java-options "{java_other_options} -Xmx{ram}" \\
  HaplotypeCaller {options} --native-pair-hmm-threads {threads} \\
  --reference {reference_sequence} \\
  --input {input} \\
  --output {output} \\
  {interval_list} {interval_padding}{intervals}{exclude_intervals}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                options=global_conf.global_get(ini_section, 'options'),
                threads=global_conf.global_get(ini_section, 'threads'),
                reference_sequence=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
                interval_list="--intervals " + str(interval_list) if interval_list else "",
                interval_padding=" \\\n --interval-padding " + str(interval_padding)  if interval_padding else "",
                input=input,
                output=output,
                intervals="".join(" \\\n  --intervals " + interval for interval in intervals) if intervals else "",
                exclude_intervals="".join(" \\\n  --exclude-intervals " + exclude_interval for exclude_interval in exclude_intervals) if exclude_intervals else ""
            )
        )

def combine_gvcf(
    inputs,
    output,
    intervals=[],
    exclude_intervals=[],
    ini_section='gatk_combine_gvcf'
    ):
    if not isinstance(inputs, list):
        inputs = [inputs]

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return gatk.combine_gvcf(inputs, output, intervals, exclude_intervals)
    else:

        return Job(
            inputs,
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
            ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CombineGVCFs {other_options} \\
  --reference {reference_sequence} \\
  {input} \\
  --output {output}{intervals}{exclude_intervals}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                other_options=global_conf.global_get(ini_section, 'other_options', required=False),
                reference_sequence=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
                input="".join(" \\\n  --variant " + input for input in inputs),
                output=output,
                intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
                exclude_intervals="".join(
                    " \\\n  --exclude-intervals " + exclude_interval for exclude_interval in exclude_intervals)
            )
        )

def GenomicsDBImport(
    inputs,
    output,
    intervals=[],
    exclude_intervals=[],
    ini_section='gatk_GenomicsDBImport'
    ):

    if not isinstance(inputs, list):
        inputs = [inputs]

    else:
        return Job(
            inputs,
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  GenomicsDBImport {other_options} \\
  {input} \\
  --output {output}{intervals}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                other_options=global_conf.global_get(ini_section, 'other_options', required=False),
                input="".join(" \\\n  --variant " + input for input in inputs),
                output=output,
                intervals="".join(" \\\n  --intervals " + interval for interval in intervals),
            )
        )

def genotype_gvcf(
    variants,
    output,
    options,
    ini_section='gatk_genotype_gvcf'
    ):

    if not isinstance(variants, list):
        variants = [variants]

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return gatk.genotype_gvcf(
            variants,
            output,
            options,
            ini_section=ini_section
        )
    else:
        return Job(
            variants,
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  GenotypeGVCFs {options} \\
  --reference {reference_sequence}{variants} \\
  --output {output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            options=options,
            reference_sequence=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
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
            interval_list=None,
            ini_section='gatk_mutect2'
            ):
    interval_padding = global_conf.global_get(ini_section, 'interval_padding')

    if interval_list:
        inputs = [inputNormal, inputTumor, interval_list]

    else:
        inputs = [inputNormal, inputTumor]

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return gatk.mutect2(inputNormal,
                            normal_name,
                            inputTumor,
                            tumor_name,
                            outputVCF,
                            interval_list)
    else:
        return Job(
            inputs,
            [outputVCF, outputVCF + ".stats", read_orientation],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
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
  --intervals {interval_list} \\
  --output {outputVCF}{pon}""".format(
        tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
        java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
        ram=global_conf.global_get(ini_section, 'ram'),
        options=global_conf.global_get(ini_section, 'options'),
        reference_sequence=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
        read_orientation=read_orientation,
        known_sites=global_conf.global_get(ini_section, 'known_sites', param_type='filepath'),
        inputNormal=inputNormal,
        normal_name=normal_name,
        inputTumor=inputTumor,
        tumor_name=tumor_name,
        outputVCF=outputVCF,
        interval_list=str(interval_list) if interval_list else "",
        interval_padding=" \\\n --interval-padding " + str(interval_padding) if interval_padding else "",
        pon=" --panel-of-normals " + global_conf.global_get(ini_section, 'pon', param_type='filepath') if global_conf.global_get(ini_section, 'pon', param_type='filepath', required=False) else ""
        )
    )

#####################
# GATK4 - Variant Filtering

def get_pileup_summaries(
    input_bam,
    output,
    ini_section='gatk_get_pileup_summaries'
    ):

    return Job(
        [input_bam],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  GetPileupSummaries {options} \\
  --input {input_bam} \\
  --variant {variants} --intervals {intervals} \\
  --output {output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            options=global_conf.global_get(ini_section, 'options'),
            variants=global_conf.global_get(ini_section, 'known_sites', param_type='filepath'),
            intervals=global_conf.global_get(ini_section, 'known_intervals', param_type='filepath'),
            input_bam=input_bam,
            output=output
        )
    )

def calculate_contamination(
        input,
        output,
        match_normal=None,
        tumor_segment=None,
        ini_section='gatk_calculate_contamination'
    ):

    return Job(
        [input, match_normal],
        [output, tumor_segment],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CalculateContamination {options} \\
  --input {input} \\
  --output {output} \\
  {normal}{tumor_segment}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            options=global_conf.global_get(ini_section, 'options'),
            input=input,
            normal=" \\\n --matched-normal " + match_normal if match_normal else "",
            tumor_segment=" \\\n --tumor-segmentation " + tumor_segment if tumor_segment else "",
            output=output,
        )
    )

def learn_read_orientation_model(
    inputs,
    output,
    ini_section='gatk_learn_read_orientation_model'
    ):

    return Job(
        inputs,
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
    command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  LearnReadOrientationModel {options} \\
  {input} \\
  --output {output}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                options=global_conf.global_get(ini_section, 'options'),
                input="".join(" \\\n  --input " + input for input in inputs),
                output=output
            )
        )

def filter_mutect_calls(
    variants,
    output,
    contamination=None,
    segment=None,
    read_orientation=None,
    ini_section='gatk_filter_mutect_calls'
    ):

    return Job(
        [variants],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  FilterMutectCalls {options} \\
  --reference {reference} \\
  --variant {variants} \\
  {contamination_table} {segment_table} {read_orientation_model} \\
  --output {output}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                options=global_conf.global_get(ini_section, 'options'),
                reference=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
                variants=variants,
                contamination_table="  \\\n --contamination-table " + contamination if contamination else "",
                segment_table="  \\\n --tumor-segmentation " + segment if segment else "",
                read_orientation_model="  \\\n --ob-priors " + read_orientation if read_orientation else "",
                output=output
            )
        )

def variant_recalibrator(
    variants,
    other_options,
    recal_output,
    tranches_output,
    R_output,
    small_sample_check=False,
    ini_section='gatk_variant_recalibrator'
    ):

    if not isinstance(variants, list):
        variants = [variants]

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return gatk.variant_recalibrator(
            variants,
            other_options,
            recal_output,
            tranches_output,
            R_output,
            small_sample_check=small_sample_check
            )
    else:
        return Job(
            variants,
            [recal_output, tranches_output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk'],
                [ini_section, 'module_R']
            ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  VariantRecalibrator {options} \\
  --reference {reference_sequence}{variants} \\
  {other_options} \\
  --output {recal_output} \\
  --tranches-file {tranches_output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            options=global_conf.global_get(ini_section, 'options'),
            reference_sequence=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
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
    apply_recal_output,
    ini_section='gatk_apply_recalibration'
    ):

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
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
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
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
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                options=global_conf.global_get(ini_section, 'options'),
                reference_sequence=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
                variants=variants,
                other_options=other_options,
                recal_input=recal_input,
                tranches_input=tranches_input,
                output=apply_recal_output
            )
        )

def variant_filtration(
    input,
    output,
    other_options,
    ini_section='gatk_variant_filtration'
    ):

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return gatk.variant_filtration(
            input,
            output,
            other_options
        )
    else:
        return Job(
            [
                input
            ],
                [
                    output
                ],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  VariantFiltration \\
  --reference {reference_sequence} \\
  --variant {variants} \\
  {other_options} \\
  --output {output}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                reference_sequence=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
                variants=input,
                other_options=other_options,
                output=output
        ))

#####################
#  Copy Number Variant Discovery
def preprocessIntervals(
    input,
    output,
    intervals,
    ini_section='gatk_processIntervals'
    ):
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk4']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  PreprocessIntervals {options} \\
  --reference {reference_sequence} \\
  --intervals {intervals} \\
  --bin-length {bin_length} \\
  --padding {padding} \\
  --output {output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            options=global_conf.global_get(ini_section, 'options'),
            reference_sequence=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
            intervals=input,
            bin_length=global_conf.global_get(ini_section, 'bin_length'),
            padding=global_conf.global_get(ini_section, 'padding'),
            output=output
        )
    )

#####################
# PICARD imported functions

def build_bam_index(
    input,
    output,
    ini_section='build_bam_index'
    ):

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return picard2.build_bam_index(
            input,
            output
        )
    else:
        return Job(
            [input],
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk4']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  BuildBamIndex \\
  VALIDATION_STRINGENCY SILENT \\
  --TMP_DIR {tmp_dir} \\
  --INPUT {input} \\
  --OUTPUT {output} """.format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                input=input,
                output=output,
            )
        )

def calculate_hs_metrics(
    input,
    output,
    intervals,
    reference_sequence=None,
    ini_section='picard_calculate_hs_metrics'
    ):

    baits_intervals = ""
    baits_intervals = global_conf.global_get(ini_section, 'baits_intervals', required=False)

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
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
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
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
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                input=input,
                output=output,
                intervals=intervals,
                baits=baits_intervals if baits_intervals != "" else intervals,
                reference_sequence=reference_sequence if reference_sequence else global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath')
            )
        )

def collect_multiple_metrics(
    input,
    output,
    reference_sequence=None,
    library_type="PAIRED_END",
    ini_section='picard_collect_multiple_metrics'
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

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
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
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk'],
                [ini_section, 'module_R']
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
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                reference_sequence=reference_sequence if reference_sequence else global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
                input=input,
                output=output,
                max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            ),
            report_files=outputs
        )

def collect_sequencing_artifacts_metrics(
    input,
    output,
    annotation_flat=None,
    reference_sequence=None,
    ini_section='picard_collect_sequencing_artifacts_metrics'
    ):

    output_dep = output + ".bait_bias_summary_metrics.txt"

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
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
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk'],
                [ini_section, 'module_R']
            ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 CollectSequencingArtifactMetrics \\
  --VALIDATION_STRINGENCY SILENT {options} \\
  --INPUT {input} \\
  --OUTPUT {output} \\
  --REFERENCE_SEQUENCE {reference} \\
  --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
            options=global_conf.global_get(ini_section, 'options'),
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            input=input,
            output=output,
            reference=reference_sequence if reference_sequence else global_conf.global_get(ini_section, 'genome_fasta'),
            max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
        )
    )

def convert_sequencing_artifacts_metrics(
    input,
    output,
    annotation_flat=None,
    reference_sequence=None,
    ini_section='picard_convert_sequencing_artifacts_metrics'
    ):

    input_dep = input + ".bait_bias_summary_metrics"

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
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
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk'],
                [ini_section, 'module_R']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  ConvertSequencingArtifactToOxoG \\
  --VALIDATION_STRINGENCY SILENT  \\
  --TMP_DIR {tmp_dir} \\
  --INPUT_BASE {input} \\
  --OUTPUT_BASE {output} \\
  --REFERENCE_SEQUENCE {reference}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                input=input,
                output=output,
                reference=reference_sequence if reference_sequence else global_conf.global_get(ini_section, 'genome_fasta'),
            )
        )

def collect_oxog_metrics(
    input,
    output,
    annotation_flat=None,
    reference_sequence=None,
    ini_section='picard_collect_oxog_metrics'
    ):

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
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
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk'],
                [ini_section, 'module_R']
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
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                input=input,
                output=output,
                dbsnp="--DB_SNP " + global_conf.global_get(ini_section, 'known_variants') if global_conf.global_get(ini_section, 'known_variants') else "",
                reference=reference_sequence if reference_sequence else global_conf.global_get(ini_section, 'genome_fasta'),
                max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            )
        )

def collect_gcbias_metrics(
    input,
    output_prefix,
    chart=None,
    summary_file=None,
    annotation_flat=None,
    reference_sequence=None,
    ini_section='picard_collect_gcbias_metrics'
    ):

    output = output_prefix +  ".gcbias_metrics.txt"
    if not chart:
        chart = output_prefix + ".gcbias_metrics.pdf"
    if not summary_file:
        summary_file = output_prefix + ".gcbias_summary_metrics.txt"
    outputs = [
        output,
        chart,
        summary_file
    ]

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
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
            outputs,
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk'],
                [ini_section, 'module_R']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CollectGcBiasMetrics \\
  --VALIDATION_STRINGENCY SILENT \\
  --ALSO_IGNORE_DUPLICATES TRUE \\
  --TMP_DIR {tmp_dir} \\
  --INPUT {input} \\
  --OUTPUT {output} \\
  --CHART_OUTPUT {chart} \\
  --SUMMARY_OUTPUT {summary_file} \\
  --REFERENCE_SEQUENCE {reference} \\
  --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            input=input,
            output=output,
            chart=chart,
            summary_file=summary_file,
            reference=reference_sequence if reference_sequence else global_conf.global_get(ini_section, 'genome_fasta'),
            max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
        ),
        report_files=outputs
    )

def fix_mate_information(
    input,
    output,
    create_index=True,
    ini_section='gatk_fix_mate_information'
    ):

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return picard2.fix_mate_information(
            input,
            output
        )
    else:
        return Job(
            [input],
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 FixMateInformation \\
 --VALIDATION_STRINGENCY SILENT \\
 {create_index} {other_options} \\
 --TMP_DIR {tmp_dir} \\
 --REFERENCE_SEQUENCE {reference} \\
 --INPUT {input} \\
 --OUTPUT {output} \\
 --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                other_options=global_conf.global_get(ini_section, 'other_options'),
                create_index="--CREATE_INDEX true" if create_index else "",
                ram=global_conf.global_get(ini_section, 'ram'), 
                reference=global_conf.global_get(ini_section, 'genome_fasta'),
                input=input,
                output=output,
                max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            ),
        )

def mark_duplicates(
    inputs,
    output,
    metrics_file,
    remove_duplicates="false",
    create_index=True,
    ini_section='gatk_mark_duplicates'
    ):

    if not isinstance(inputs, list):
        inputs = [inputs]

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return picard2.mark_duplicates(
            inputs,
            output,
            metrics_file,
            remove_duplicates,
            create_index,
            ini_section=ini_section
        )
    else:
        outputs = [output, metrics_file]
        if create_index:
            outputs.append(re.sub(r'\b(bam|cram)\b',
                                  lambda match: match.group() + (".bai" if match.group() == "bam" else ".crai"),
                                  output)
                           )
        return Job(
            inputs,
            outputs,
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
            ],
            command="""\
rm -rf {output}.part && \\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 MarkDuplicates {other_options} \\
 --REFERENCE_SEQUENCE {reference} \\
 --REMOVE_DUPLICATES {remove_duplicates} \\
 --VALIDATION_STRINGENCY SILENT \\
 {create_index} \\
 --TMP_DIR {tmp_dir} \\
 {inputs} \\
 --OUTPUT {output} \\
 --METRICS_FILE {metrics_file} \\
 --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                other_options=global_conf.global_get(ini_section, 'other_options', required=False),
                reference=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                remove_duplicates=remove_duplicates,
                create_index="--CREATE_INDEX true" if create_index else "",
                inputs=" \\\n  ".join("--INPUT " + input for input in inputs),
                output=output,
                metrics_file=metrics_file,
                max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            ),
        )

def mark_duplicates_mate_cigar(
    inputs,
    output,
    metrics_file,
    remove_duplicates="false",
    ini_section='mark_duplicates_mate_cigar'
    ):

    if not isinstance(inputs, list):
        inputs = [inputs]
    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return picard2.mark_duplicates_mate_cigar(
            inputs,
            output,
            metrics_file,
            remove_duplicates
        )
    else:
        return Job(
            inputs,
            [output, re.sub(r"\.([sb])am$", ".\\1ai", output), metrics_file],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
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
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                remove_duplicates=remove_duplicates,
                inputs=" \\\n  ".join("--INPUT " + input for input in inputs),
                output=output,
                metrics_file=metrics_file,
                max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            )
        )

def picard_mark_duplicates_mate_cigar(
    inputs,
    output,
    metrics_file,
    remove_duplicates="false",
    ini_section='picard_mark_duplicates_mate_cigar'
    ):

    if not isinstance(inputs, list):
        inputs = [inputs]
    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return picard2.mark_duplicates_mate_cigar(
            inputs,
            output,
            metrics_file,
            remove_duplicates
        )
    else:
        return Job(
            inputs,
            [output, re.sub(r"\.([sb])am$", ".\\1ai", output), metrics_file],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
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
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                remove_duplicates=remove_duplicates,
                inputs=" \\\n  ".join("--INPUT " + input for input in inputs),
                output=output,
                metrics_file=metrics_file,
                max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            ),
        )

def merge_sam_files(
    inputs,
    output,
    ini_section='picard_merge_sam_files'
    ):

    if not isinstance(inputs, list):
        inputs = [inputs]
    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return picard2.merge_sam_files(
            inputs,
            output
        )
    else:
        return Job(
            inputs,
            [
                output,
                re.sub(r"\.([sb])am$", ".\\1ai", output)
            ],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk']
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
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                inputs=" \\\n ".join(["--INPUT " + input for input in inputs]),
                output=output,
                max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            ),
        )

# Reorder BAM/SAM files based on reference/dictionary
def reorder_sam(
    input,
    output,
    ini_section='picard_reorder_sam'
    ):

    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
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
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                input=input,
                output=output,
                reference=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
                max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            ),
        )

# Convert SAM/BAM file to fastq format
def sam_to_fastq(
    input,
    fastq,
    second_end_fastq=None,
    ini_section='gatk_sam_to_fastq'
    ):

    return Job(
        [input],
        [
            fastq,
            second_end_fastq
        ],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
 SamToFastq {other_options} \\
 --CREATE_MD5_FILE true \\
 --INPUT {input} \\
 --FASTQ {fastq}{second_end_fastq}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                other_options=global_conf.global_get(ini_section, 'other_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
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

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
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
                re.sub(r"\.([sb])am$", ".\\1ai", output) if sort_order == "coordinate" else None
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
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                input=input,
                output=output,
                sort_order=sort_order,
                max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            ),
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
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                inputs=" \\\n  ".join(["INPUT " + input for input in inputs]),
                output=output,
                seq_dict=global_conf.global_get(ini_section, 'genome_dictionary', param_type='filepath')
            )
        )

def collect_rna_metrics(
    input,
    output,
    annotation_flat=None,
    reference_sequence=None,
    ini_section='picard_collect_rna_metrics'
    ):
    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return picard2.collect_rna_metrics(
            input,
            output,
            annotation_flat=None,
            reference_sequence=None
        )
    else:
        return Job(
            [input],
            # collect specific RNA metrics (exon rate, strand specificity, etc...)
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk'],
                [ini_section, 'module_R']
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
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                input=input,
                output=output,
                ref_flat=annotation_flat if annotation_flat else global_conf.global_get(ini_section, 'annotation_flat'),
                strand_specificity=global_conf.global_get(ini_section, 'strand_info'),
                min_length=global_conf.global_get(ini_section, 'minimum_length', param_type='int'),
                reference=reference_sequence if reference_sequence else global_conf.global_get(ini_section, 'genome_fasta'),
                max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            )
        )


def crosscheck_fingerprint(
    inputs,
    output,
    ini_section='gatk_crosscheck_fingerprint'
    ):

    matrix=output + ".matrix"
    return Job(
        inputs,
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CrosscheckFingerprints {options} \\
  --REFERENCE_SEQUENCE {reference} \\
  --VALIDATION_STRINGENCY SILENT \\
  --TMP_DIR {tmp_dir} \\
  {inputs} \\
  --HAPLOTYPE_MAP {haplotype_database} \\
  --LOD_THRESHOLD {lod_threshold} \\
  --OUTPUT {output} \\
  --MATRIX_OUTPUT {matrix}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                options=global_conf.global_get(ini_section, 'options'),
                reference = global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                haplotype_database=global_conf.global_get(ini_section, 'haplotype_database'),
                lod_threshold=global_conf.global_get(ini_section, 'lod_threshold'),
                ram=global_conf.global_get(ini_section, 'ram'),
                inputs=" \\\n  ".join("--INPUT " + str(input) for input in inputs),
                output=output,
                matrix=matrix,
            )
        )

def cluster_crosscheck_metrics(
    input,
    output,
    ini_section='gatk_cluster_crosscheck_metrics'
    ):
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  ClusterCrosscheckMetrics {options} \\
  --VALIDATION_STRINGENCY SILENT \\
  --TMP_DIR {tmp_dir} \\
  --INPUT {input} \\
  --LOD_THRESHOLD {lod_threshold} \\
  --OUTPUT {output} """.format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            options=global_conf.global_get(ini_section, 'options'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            lod_threshold=global_conf.global_get(ini_section, 'lod_threshold'),
            ram=global_conf.global_get(ini_section, 'ram'),
            input=input,
            output=output,
        )
    )
def scatterIntervalsByNs(
    reference,
    output,
    ini_section='gatk_scatterIntervalsByNs'
    ):

    return Job(
        [reference],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  ScatterIntervalsByNs \\
  {options} \\
  --OUTPUT_TYPE {output_type} \\
  --REFERENCE {reference} \\
  --OUTPUT {output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            options=global_conf.global_get(ini_section, 'options'),
            output_type=global_conf.global_get(ini_section, 'output_type'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            reference=reference,
            output=output
        )
    )

def bed2interval_list(
    dictionary,
    bed,
    output,
    ini_section='gatk_bed2interval_list'
    ):

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return gatk.bed2interval_list(
            dictionary,
            bed,
            output
            )

    return Job(
        [dictionary, bed],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  BedToIntervalList \\
  --INPUT {bed} \\
  --SEQUENCE_DICTIONARY {dictionary} \\
  --OUTPUT {output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            dictionary=dictionary if dictionary else global_conf.global_get(ini_section, 'genome_dictionary', param_type='filepath'),
            bed=bed,
            output=output
            )
        )

def interval_list2bed(
    input,
    output,
    ini_section='gatk_interval_list2bed'
    ):
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  IntervalListToBed \\
  --INPUT {input} \\
  --OUTPUT {output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            input=input,
            output=output
            )
        )

def preProcessInterval(
    reference,
    intervals,
    output,
    options=None,
    ini_section='gatk_preProcessInterval'
    ):

    return Job(
        [intervals],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  PreprocessIntervals {options} \\
  --reference {reference} \\
  --intervals {intervals} \\
  --output {output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            options=options if options else global_conf.global_get(ini_section, 'options'),
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            reference=reference if reference else global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
            intervals=intervals,
            output=output
        )
    )

def splitInterval(
    intervals,
    output,
    jobs,
    options=None,
    ini_section='gatk_splitInterval'
    ):

    interval_list = []
    for idx in range(jobs):
        interval_list.append(os.path.join(output, str(idx).zfill(4) + "-scattered.interval_list"))

    return Job(
        [intervals],
        interval_list,
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_gatk']
        ],
        command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  SplitIntervals {options} \\
  --scatter-count {jobs} \\
  --reference {reference} \\
  --intervals {intervals} \\
  --output {output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            options=options if options else global_conf.global_get(ini_section, 'options'),
            jobs=jobs,
            java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            reference=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
            intervals=intervals,
            output=output
        )
    )

def collect_wgs_metrics(
    input,
    output,
    reference_sequence=None,
    ini_section='picard_collect_wgs_metrics'
    ):

    if global_conf.global_get(ini_section, 'module_gatk').split("/")[2] < "4":
        return picard2.collect_wgs_metrics(
            input,
            output,
            reference_sequence
        )
    else:
        return Job(
            [input],
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_gatk'],
            ],
            command="""\
gatk --java-options "-Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram}" \\
  CollectWgsMetrics {options} \\
  --VALIDATION_STRINGENCY SILENT  \\
  --TMP_DIR {tmp_dir} \\
  --INPUT {input} \\
  --OUTPUT {output} \\
  --REFERENCE_SEQUENCE {reference} \\
  --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'gatk_java_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                options=global_conf.global_get(ini_section, 'options'),
                input=input,
                output=output,
                reference=reference_sequence if reference_sequence else global_conf.global_get(ini_section, 'genome_fasta'),
                max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            )
        )

def parse_bases_over_q30_percent_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export bases_over_q30_percent=`awk 'BEGIN {{FS="\t"}}; {{if ($1 ~ /^[0-9]+/) {{if ($1<30) {{below+=$2}} else if ($1>=30) {{above+=$2}}}}}} END {{printf "%.0f", 100*above/(above+below)}}' {input_file}`"""
        )
def parse_duplicate_rate_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export duplication_percent=`grep -vE "^$" {input_file} | awk '/^## METRICS CLASS/{{flag=1; next}} /## HISTOGRAM/{{flag=0}} flag && NF>1 {{sum+=$9; count++}} END {{if(count>0) print sum/(count-1)}}'`"""
        )
def parse_mean_insert_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export mean_insert_size=`grep -A1 "^MEDIAN" {input_file} | cut -f6 | grep -vE "MEAN|^$"`"""
        )

def parse_stdev_insert_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export stdev_insert_size=`grep -A1 "^MEDIAN" {input_file} | cut -f7 | grep -vE "STANDARD|^$"`"""
        )

def parse_mode_insert_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export mode_insert_size=`grep -A1 "^MEDIAN" {input_file} | cut -f2 | grep -vE "MODE|^$"`"""
        )

def parse_total_read_pairs_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export total_read_pairs=`grep "^PAIR" {input_file} | cut -f2`"""
        )

def parse_aligned_pairs_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export aligned_pairs_percent=`grep "^PAIR" {input_file} | cut -f7`"""
        )

def parse_high_quality_read_pairs_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export hq_read_pairs=`grep "^PAIR" {input_file} | cut -f9`"""
        )

def parse_chimeras_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export chimeras_percent=`grep "^PAIR" {input_file} | cut -f28`"""
        )

def parse_bed_bait_set_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export bed_bait_set=`grep -A1 "^BAIT" {input_file} | cut -f1 | grep -vE "BAIT|^$"`"""
        )

def parse_off_target_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export off_target_percent=`grep -A1 "^BAIT" {input_file} | cut -f8 | grep -vE "PCT|^$"`"""
        )

def parse_target_basepair_size_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export target_basepair_size=`grep -A1 "^BAIT" {input_file} | cut -f21 | grep -vE "TARGET|^$"`"""
        )

def parse_total_reads_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export total_reads=`grep -A1 "^BAIT" {input_file} | cut -f23 | grep -vE "TOTAL|^$"`"""
        )

def parse_dedup_reads_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export dedup_reads=`grep -A1 "^BAIT" {input_file} | cut -f26 | grep -vE "PF|^$"`"""
        )

def parse_mean_target_coverage_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export mean_target_coverage=`grep -A1 "^BAIT" {input_file} | cut -f34 | grep -vE "MEAN|^$"`"""
        )

def parse_median_target_coverage_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export median_target_coverage=`grep -A1 "^BAIT" {input_file} | cut -f35 | grep -vE "MEDIAN|^$"`"""
        )

def parse_dup_rate_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export duplicate_rate_percent=`grep -A1 "^BAIT" {input_file} | cut -f39 | grep -vE "PCT|^$"`"""
        )

def parse_low_mapping_rate_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export low_mapping_rate_percent=`grep -A1 "^BAIT" {input_file} | cut -f41 | grep -vE "PCT|^$"`"""
        )

def parse_read_overlap_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export read_overlap_percent=`grep -A1 "^BAIT" {input_file} | cut -f43 | grep -vE "PCT|^$"`"""
        )