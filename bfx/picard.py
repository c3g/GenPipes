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

def build_bam_index(input, output):

    return Job(
        [input],
        [output],
        [
            ['build_bam_index', 'module_java'],
            ['build_bam_index', 'module_picard']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/BuildBamIndex.jar \\
  VALIDATION_STRINGENCY=SILENT \\
  INPUT={input} \\
  OUTPUT={output} """.format(
        tmp_dir=config.param('build_bam_index', 'tmp_dir'),
        java_other_options=config.param('build_bam_index', 'java_other_options'),
        ram=config.param('build_bam_index', 'ram'),
        input=input,
        output=output,
        )
    )

def calculate_hs_metrics(input, output, intervals, reference_sequence=None):
    
    baits_intervals = ""
    baits_intervals = config.param('picard_calculate_hs_metrics', 'baits_intervals', required = False)

    return Job(
        [input, intervals],
        [output],
        [
            ['picard_calculate_hs_metrics', 'module_java'],
            ['picard_calculate_hs_metrics', 'module_picard']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/CalculateHsMetrics.jar \\
  TMP_DIR={tmp_dir} \\
  INPUT={input} \\
  OUTPUT={output} \\
  BAIT_INTERVALS={baits} \\
  TARGET_INTERVALS={intervals} \\
  REFERENCE_SEQUENCE={reference_sequence}""".format(
        tmp_dir=config.param('picard_calculate_hs_metrics', 'tmp_dir'),
        java_other_options=config.param('picard_calculate_hs_metrics', 'java_other_options'),
        ram=config.param('picard_calculate_hs_metrics', 'ram'),
        input=input,
        output=output,
        intervals=intervals,
        baits=baits_intervals if baits_intervals != "" else intervals,
        reference_sequence=reference_sequence if reference_sequence else config.param('picard_calculate_hs_metrics', 'genome_fasta', type='filepath')
        )
    )

def collect_multiple_metrics(input, output, reference_sequence=None , library_type="PAIRED_END"):
    
    if  library_type == "PAIRED_END" :
        outputs = [
         output + ".quality_by_cycle.pdf",
         output + ".alignment_summary_metrics",
         output + ".insert_size_Histogram.pdf",
         output + ".insert_size_metrics",
         output + ".quality_by_cycle_metrics",
         output + ".quality_distribution_metrics",
         output + ".quality_distribution.pdf"
        ]
    else :
        outputs = [
         output + ".quality_by_cycle.pdf",
         output + ".alignment_summary_metrics",
         output + ".quality_by_cycle_metrics",
         output + ".quality_distribution_metrics",
         output + ".quality_distribution.pdf"
        ]

    return Job(
        [input],
        outputs,
        [
            ['picard_collect_multiple_metrics', 'module_java'],
            ['picard_collect_multiple_metrics', 'module_picard'],
            ['picard_collect_multiple_metrics', 'module_R']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/CollectMultipleMetrics.jar \\
  PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \\
  TMP_DIR={tmp_dir} \\
  REFERENCE_SEQUENCE={reference_sequence} \\
  INPUT={input} \\
  OUTPUT={output} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=config.param('picard_collect_multiple_metrics', 'tmp_dir'),
        java_other_options=config.param('picard_collect_multiple_metrics', 'java_other_options'),
        ram=config.param('picard_collect_multiple_metrics', 'ram'),
        reference_sequence=reference_sequence if reference_sequence else config.param('picard_collect_multiple_metrics', 'genome_fasta', type='filepath'),
        input=input,
        output=output,
        max_records_in_ram=config.param('picard_collect_multiple_metrics', 'max_records_in_ram', type='int')
        )
    )

def fix_mate_information(input, output):

    return Job(
        [input],
        [output],
        [
            ['fixmate', 'module_java'],
            ['fixmate', 'module_picard']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/FixMateInformation.jar \\
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate \\
  TMP_DIR={tmp_dir} \\
  INPUT={input} \\
  OUTPUT={output} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=config.param('picard_fix_mate_information', 'tmp_dir'),
        java_other_options=config.param('picard_fix_mate_information', 'java_other_options'),
        ram=config.param('picard_fix_mate_information', 'ram'),
        input=input,
        output=output,
        max_records_in_ram=config.param('picard_fix_mate_information', 'max_records_in_ram', type='int')
        ),
        removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
    )

def mark_duplicates(inputs, output, metrics_file):

    return Job(
        inputs,
        [output, re.sub("\.([sb])am$", ".\\1ai", output), metrics_file],
        [
            ['picard_mark_duplicates', 'module_java'],
            ['picard_mark_duplicates', 'module_picard']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/MarkDuplicates.jar \\
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
  TMP_DIR={tmp_dir} \\
  {inputs} \\
  OUTPUT={output} \\
  METRICS_FILE={metrics_file} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=config.param('picard_mark_duplicates', 'tmp_dir'),
        java_other_options=config.param('picard_mark_duplicates', 'java_other_options'),
        ram=config.param('picard_mark_duplicates', 'ram'),
        inputs=" \\\n  ".join(["INPUT=" + input for input in inputs]),
        output=output,
        metrics_file=metrics_file,
        max_records_in_ram=config.param('picard_mark_duplicates', 'max_records_in_ram', type='int')
        ),
        removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
    )

def merge_sam_files(inputs, output):

    return Job(
        inputs,
        [output, re.sub("\.([sb])am$", ".\\1ai", output)],
        [
            ['picard_merge_sam_files', 'module_java'],
            ['picard_merge_sam_files', 'module_picard']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/MergeSamFiles.jar \\
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \\
  TMP_DIR={tmp_dir} \\
  {inputs} \\
  OUTPUT={output} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=config.param('picard_merge_sam_files', 'tmp_dir'),
        java_other_options=config.param('picard_merge_sam_files', 'java_other_options'),
        ram=config.param('picard_merge_sam_files', 'ram'),
        inputs=" \\\n  ".join(["INPUT=" + input for input in inputs]),
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
            ['reorder_sam', 'module_picard']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/ReorderSam.jar \\
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
  TMP_DIR={tmp_dir} \\
  INPUT={input} \\
  OUTPUT={output} \\
  REFERENCE={reference} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
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
            ['picard_sam_to_fastq', 'module_picard']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/SamToFastq.jar \\
  VALIDATION_STRINGENCY=LENIENT \\
  INPUT={input} \\
  FASTQ={fastq}{second_end_fastq}""".format(
        tmp_dir=config.param('picard_sam_to_fastq', 'tmp_dir'),
        java_other_options=config.param('picard_sam_to_fastq', 'java_other_options'),
        ram=config.param('picard_sam_to_fastq', 'ram'),
        input=input,
        fastq=fastq,
        second_end_fastq=" \\\n  SECOND_END_FASTQ=" + second_end_fastq if second_end_fastq else ""
        ),
        removable_files=[fastq, second_end_fastq]
    )

def sort_sam(input, output, sort_order="coordinate",ini_section='picard_sort_sam'):

    return Job(
        [input],
        # Add SAM/BAM index as output only when writing a coordinate-sorted BAM file
        [output, re.sub("\.([sb])am$", ".\\1ai", output) if sort_order == "coordinate" else None],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_picard']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/SortSam.jar \\
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
  TMP_DIR={tmp_dir} \\
  INPUT={input} \\
  OUTPUT={output} \\
  SORT_ORDER={sort_order} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
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

    return Job(
        inputs,
        # Add SAM/BAM index as output only when writing a coordinate-sorted BAM file
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_picard']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/SortVcf.jar \\
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

def collect_rna_metrics(input, output, annotation_flat=None,reference_sequence=None):

    return Job(
        [input],
        # collect specific RNA metrics (exon rate, strand specificity, etc...)
        [output],
        [
            ['picard_collect_rna_metrics', 'module_java'],
            ['picard_collect_rna_metrics', 'module_picard'],
            ['picard_collect_rna_metrics', 'module_R']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/CollectRnaSeqMetrics.jar \\
  VALIDATION_STRINGENCY=SILENT  \\
  TMP_DIR={tmp_dir} \\
  INPUT={input} \\
  OUTPUT={output} \\
  REF_FLAT={ref_flat} \\
  STRAND_SPECIFICITY={strand_specificity} \\
  MINIMUM_LENGTH={min_length} \\
  REFERENCE_SEQUENCE={reference} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=config.param('picard_collect_rna_metrics', 'tmp_dir'),
        java_other_options=config.param('picard_collect_rna_metrics', 'java_other_options'),
        ram=config.param('picard_collect_rna_metrics', 'ram'),
        input=input,
        output=output,
        ref_flat=annotation_flat if annotation_flat else config.param('picard_collect_rna_metrics', 'annotation_flat'),
        strand_specificity=config.param('picard_collect_rna_metrics', 'strand_info'),
        min_length=config.param('picard_collect_rna_metrics', 'minimum_length',type='int'),
        reference=reference_sequence if reference_sequence else config.param('picard_collect_rna_metrics', 'genome_fasta'),
        max_records_in_ram=config.param('picard_collect_rna_metrics', 'max_records_in_ram', type='int')
        )
    )
