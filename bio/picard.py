#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def calculate_hs_metrics(input, output, intervals, reference_sequence=None):

    job = Job([input], [output], [['calculate_hs_metrics', 'moduleVersion.java'], ['calculate_hs_metrics', 'moduleVersion.picard']])

    if not reference_sequence:
        reference_sequence = config.param('calculate_hs_metrics', 'referenceFasta', type='filepath')

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$PICARD_HOME/CalculateHsMetrics.jar \\
  TMP_DIR={tmp_dir} \\
  INPUT={input} \\
  OUTPUT={output} \\
  BAIT_INTERVALS={intervals} \\
  TARGET_INTERVALS={intervals} \\
  REFERENCE_SEQUENCE={reference_sequence}""".format(
        tmp_dir=config.param('calculate_hs_metrics', 'tmpDir'),
        extra_java_flags=config.param('calculate_hs_metrics', 'extraJavaFlags'),
        ram=config.param('calculate_hs_metrics', 'ram'),
        input=input,
        output=output,
        intervals=intervals,
        reference_sequence=reference_sequence
    )

    return job

def collect_multiple_metrics(input, output, reference_sequence=None):

    job = Job(
        [input],
        [output + ".quality_by_cycle.pdf"],
        [
            ['collectMetrics', 'moduleVersion.java'],
            ['collectMetrics', 'moduleVersion.picard'],
            ['collectMetrics', 'moduleVersion.R']
        ]
    )

    if not reference_sequence:
        reference_sequence = config.param('collectMetrics', 'referenceFasta', type='filepath')

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$PICARD_HOME/CollectMultipleMetrics.jar \\
  PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \\
  TMP_DIR={tmp_dir} \\
  REFERENCE_SEQUENCE={reference_sequence} \\
  INPUT={input} \\
  OUTPUT={output} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=config.param('collectMetrics', 'tmpDir'),
        extra_java_flags=config.param('collectMetrics', 'extraJavaFlags'),
        ram=config.param('collectMetrics', 'collectMetricsRam'),
        reference_sequence=reference_sequence,
        input=input,
        output=output,
        max_records_in_ram=config.param('collectMetrics', 'collectMetricsRecInRam', type='int')
    )

    return job

def fix_mate_information(input, output):

    job = Job([input], [output], [['fixmate', 'moduleVersion.java'], ['fixmate', 'moduleVersion.picard']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$PICARD_HOME/FixMateInformation.jar \\
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate \\
  TMP_DIR={tmp_dir} \\
  INPUT={input} \\
  OUTPUT={output} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=config.param('fixmate', 'tmpDir'),
        extra_java_flags=config.param('fixmate', 'extraJavaFlags'),
        ram=config.param('fixmate', 'fixmateRam'),
        input=input,
        output=output,
        max_records_in_ram=config.param('fixmate', 'fixmateRecInRam', type='int')
    )

    return job

def mark_duplicates(inputs, output, metrics_file):

    job = Job(inputs, [output, metrics_file], [['markDup', 'moduleVersion.java'], ['markDup', 'moduleVersion.picard']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$PICARD_HOME/MarkDuplicates.jar \\
  REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
  TMP_DIR={tmp_dir} \\
  {inputs} \\
  OUTPUT={output} \\
  METRICS_FILE={metrics_file} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=config.param('markDup', 'tmpDir'),
        extra_java_flags=config.param('markDup', 'extraJavaFlags'),
        ram=config.param('markDup', 'markDupRam'),
        inputs=" \\\n  ".join(["INPUT=" + input for input in inputs]),
        output=output,
        metrics_file=metrics_file,
        max_records_in_ram=config.param('markDup', 'markDupRecInRam', type='int')
    )

    return job

def merge_sam_files(inputs, output):

    job = Job(inputs, [output], [['mergeSamFiles', 'moduleVersion.java'], ['mergeSamFiles', 'moduleVersion.picard']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$PICARD_HOME/MergeSamFiles.jar \\
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \\
  TMP_DIR={tmp_dir} \\
  {inputs} \\
  OUTPUT={output} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=config.param('mergeSamFiles', 'tmpDir'),
        extra_java_flags=config.param('mergeSamFiles', 'extraJavaFlags'),
        ram=config.param('mergeSamFiles', 'mergeRam'),
        inputs=" \\\n  ".join(["INPUT=" + input for input in inputs]),
        output=output,
        max_records_in_ram=config.param('mergeSamFiles', 'mergeRecInRam', type='int')
    )

    return job

# Reorder BAM/SAM files based on reference/dictionary
def reorder_sam(input, output):

    job = Job([input], [output], [['reorderSam', 'moduleVersion.java'], ['reorderSam', 'moduleVersion.picard']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$PICARD_HOME/ReorderSam.jar \\
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
  TMP_DIR={tmp_dir} \\
  INPUT={input} \\
  OUTPUT={output} \\
  REFERENCE={reference} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=config.param('reorderSam', 'tmpDir'),
        extra_java_flags=config.param('reorderSam', 'extraJavaFlags'),
        ram=config.param('reorderSam', 'reorderRam'),
        input=input,
        output=output,
        reference=config.param('reorderSam', 'referenceFasta', type='filepath'),
        max_records_in_ram=config.param('reorderSam', 'reorderRecInRam', type='int')
    )

    return job

# Convert SAM/BAM file to fastq format
def sam_to_fastq(input, fastq, second_end_fastq):

    job = Job(
        [input],
        [fastq, second_end_fastq],
        [
            ['samToFastq', 'moduleVersion.java'],
            ['samToFastq', 'moduleVersion.picard']
        ]
    )

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$PICARD_HOME/SamToFastq.jar \\
  INPUT={input} \\
  FASTQ={fastq}{second_end_fastq}""".format(
        tmp_dir=config.param('samToFastq', 'tmpDir'),
        extra_java_flags=config.param('samToFastq', 'extraJavaFlags'),
        ram=config.param('samToFastq', 'samToFastqRam'),
        input=input,
        fastq=fastq,
        second_end_fastq=" \\\n  SECOND_END_FASTQ=" + second_end_fastq if second_end_fastq else ""
    )

    return job

def sort_sam(input, output, sort_order):

    job = Job([input], [output], [['sortSam', 'moduleVersion.java'], ['sortSam', 'moduleVersion.picard']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$PICARD_HOME/SortSam.jar \\
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
  TMP_DIR={tmp_dir} \\
  INPUT={input} \\
  OUTPUT={output} \\
  SORT_ORDER={sort_order} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=config.param('sortSam', 'tmpDir'),
        extra_java_flags=config.param('sortSam', 'extraJavaFlags'),
        ram=config.param('sortSam', 'sortRam'),
        input=input,
        output=output,
        sort_order = sort_order,
        max_records_in_ram=config.param('sortSam', 'sortRecInRam', type='int')
    )

    return job
