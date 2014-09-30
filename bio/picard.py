#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def calculate_hs_metrics(input, output, intervals, reference_sequence=None):

    job = Job([input], [output], [['picard_calculate_hs_metrics', 'module_java'], ['picard_calculate_hs_metrics', 'module_picard']])

    job.command = """\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/CalculateHsMetrics.jar \\
  TMP_DIR={tmp_dir} \\
  INPUT={input} \\
  OUTPUT={output} \\
  BAIT_INTERVALS={intervals} \\
  TARGET_INTERVALS={intervals} \\
  REFERENCE_SEQUENCE={reference_sequence}""".format(
        tmp_dir=config.param('picard_calculate_hs_metrics', 'tmp_dir'),
        java_other_options=config.param('picard_calculate_hs_metrics', 'java_other_options'),
        ram=config.param('picard_calculate_hs_metrics', 'ram'),
        input=input,
        output=output,
        intervals=intervals,
        reference_sequence=reference_sequence if reference_sequence else config.param('picard_calculate_hs_metrics', 'genome_fasta', type='filepath')
    )

    return job

def collect_multiple_metrics(input, output, reference_sequence=None):

    job = Job(
        [input],
        [
         output + ".quality_by_cycle.pdf",
         output + ".alignment_summary_metrics",
         output + ".insert_size_histogram.pdf",
         output + ".insert_size_metrics",
         output + ".quality_by_cycle_metrics",
         output + ".quality_distribution_metrics",
         output + ".quality_distribution.pdf"
        ],
        [
            ['picard_collect_multiple_metrics', 'module_java'],
            ['picard_collect_multiple_metrics', 'module_picard'],
            ['picard_collect_multiple_metrics', 'module_R']
        ]
    )

    job.command = """\
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

    return job

def fix_mate_information(input, output):

    job = Job([input], [output], [['fixmate', 'module_java'], ['fixmate', 'module_picard']])

    job.command = """\
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
    )

    return job

def mark_duplicates(inputs, output, metrics_file):

    job = Job(inputs, [output, metrics_file], [['picard_mark_duplicates', 'module_java'], ['picard_mark_duplicates', 'module_picard']])

    job.command = """\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/MarkDuplicates.jar \\
  REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
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
    )

    return job

def merge_sam_files(inputs, output):

    job = Job(inputs, [output], [['picard_merge_sam_files', 'module_java'], ['picard_merge_sam_files', 'module_picard']])

    job.command = """\
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
    )

    return job

# Reorder BAM/SAM files based on reference/dictionary
def reorder_sam(input, output):

    job = Job([input], [output], [['reorder_sam', 'module_java'], ['reorder_sam', 'module_picard']])

    job.command = """\
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
    )

    return job

# Convert SAM/BAM file to fastq format
def sam_to_fastq(input, fastq, second_end_fastq):

    job = Job(
        [input],
        [fastq, second_end_fastq],
        [
            ['picard_sam_to_fastq', 'module_java'],
            ['picard_sam_to_fastq', 'module_picard']
        ]
    )

    job.command = """\
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
    )

    return job

def sort_sam(input, output, sort_order):

    job = Job([input], [output], [['picard_sort_sam', 'module_java'], ['picard_sort_sam', 'module_picard']])

    job.command = """\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/SortSam.jar \\
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
  TMP_DIR={tmp_dir} \\
  INPUT={input} \\
  OUTPUT={output} \\
  SORT_ORDER={sort_order} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=config.param('picard_sort_sam', 'tmp_dir'),
        java_other_options=config.param('picard_sort_sam', 'java_other_options'),
        ram=config.param('picard_sort_sam', 'ram'),
        input=input,
        output=output,
        sort_order=sort_order,
        max_records_in_ram=config.param('picard_sort_sam', 'max_records_in_ram', type='int')
    )

    return job
