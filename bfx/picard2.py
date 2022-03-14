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

# MUGQIC Modules
from core.config import global_config_parser
from core.job import Job
from . import picard
from . import gatk4

def build_bam_index(input, output):

    if global_config_parser.param('build_bam_index', 'module_picard').split("/")[2] < "2":
        return picard.build_bam_index(input, output)
    else:
        return Job(
            [input],
            [output],
            [
                ['build_bam_index', 'module_java'],
                ['build_bam_index', 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar BuildBamIndex \\
 VALIDATION_STRINGENCY=SILENT \\
 INPUT={input} \\
 OUTPUT={output} """.format(
            tmp_dir=global_config_parser.param('build_bam_index', 'tmp_dir'),
            java_other_options=global_config_parser.param('build_bam_index', 'java_other_options'),
            ram=global_config_parser.param('build_bam_index', 'ram'),
            input=input,
            output=output,
            )
        )

def calculate_hs_metrics(input, output, intervals, reference_sequence=None):

    baits_intervals = ""
    baits_intervals = global_config_parser.param('picard_calculate_hs_metrics', 'baits_intervals', required = False)

    if global_config_parser.param('picard_calculate_hs_metrics', 'module_picard').split("/")[2] < "2":
        return picard.calculate_hs_metrics(input, output, intervals, reference_sequence)
    else:
        return Job(
            [input, intervals],
            [output],
            [
                ['picard_calculate_hs_metrics', 'module_java'],
                ['picard_calculate_hs_metrics', 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar CollectHsMetrics \\
 TMP_DIR={tmp_dir} \\
 INPUT={input} \\
 OUTPUT={output} \\
 BAIT_INTERVALS={baits} \\
 TARGET_INTERVALS={intervals} \\
 REFERENCE_SEQUENCE={reference_sequence}""".format(
            tmp_dir=global_config_parser.param('picard_calculate_hs_metrics', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_calculate_hs_metrics', 'java_other_options'),
            ram=global_config_parser.param('picard_calculate_hs_metrics', 'ram'),
            input=input,
            output=output,
            intervals=intervals,
            baits=baits_intervals if baits_intervals != "" else intervals,
            reference_sequence=reference_sequence if reference_sequence else global_config_parser.param('picard_calculate_hs_metrics', 'genome_fasta', param_type='filepath')
            )
        )

def collect_multiple_metrics(input, output, reference_sequence=None, library_type="PAIRED_END"):

    if  library_type == "PAIRED_END" :
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
    else :
        outputs = [
         output + ".quality_by_cycle.pdf",
         output + ".alignment_summary_metrics",
         output + ".quality_by_cycle_metrics",
         output + ".quality_distribution_metrics",
         output + ".quality_distribution.pdf"
        ]

    if global_config_parser.param('picard_collect_multiple_metrics', 'module_picard').split("/")[2] < "2":
        return picard.collect_multiple_metrics(input, output, reference_sequence, library_type)
    
    else:
        
        return Job(
            [input],
            outputs,
            [
                ['picard_collect_multiple_metrics', 'module_java'],
                ['picard_collect_multiple_metrics', 'module_picard'],
                ['picard_collect_multiple_metrics', 'module_R']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \\
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \\
 TMP_DIR={tmp_dir} \\
 REFERENCE_SEQUENCE={reference_sequence} \\
 INPUT={input} \\
 OUTPUT={output} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
            tmp_dir=global_config_parser.param('picard_collect_multiple_metrics', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_collect_multiple_metrics', 'java_other_options'),
            ram=global_config_parser.param('picard_collect_multiple_metrics', 'ram'),
            reference_sequence=reference_sequence if reference_sequence else global_config_parser.param('picard_collect_multiple_metrics', 'genome_fasta', param_type='filepath'),
            input=input,
            output=output,
            max_records_in_ram=global_config_parser.param('picard_collect_multiple_metrics', 'max_records_in_ram', param_type='int')
            )
        )
def collect_sequencing_artifacts_metrics(input, output, annotation_flat=None,reference_sequence=None):
        output_dep = output + ".bait_bias_summary_metrics"

        return Job(
            [input],
            [output_dep],
            [
                ['picard_collect_sequencing_artifacts_metrics', 'module_java'],
                ['picard_collect_sequencing_artifacts_metrics', 'module_picard'],
                ['picard_collect_sequencing_artifacts_metrics', 'module_R']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar CollectSequencingArtifactMetrics \\
 VALIDATION_STRINGENCY=SILENT {options} \\
 TMP_DIR={tmp_dir} \\
 INPUT={input} \\
 OUTPUT={output} \\
 REFERENCE_SEQUENCE={reference} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
            options=global_config_parser.param('picard_collect_sequencing_artifacts_metrics', 'options'),
            tmp_dir=global_config_parser.param('picard_collect_sequencing_artifacts_metrics', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_collect_sequencing_artifacts_metrics', 'java_other_options'),
            ram=global_config_parser.param('picard_collect_sequencing_artifacts_metrics', 'ram'),
            input=input,
            output=output,
            reference=reference_sequence if reference_sequence else global_config_parser.param('picard_collect_sequencing_artifacts_metrics', 'genome_fasta'),
            max_records_in_ram=global_config_parser.param('picard_collect_sequencing_artifacts_metrics', 'max_records_in_ram', param_type='int')
            )
        )

def convert_sequencing_artifacts_metrics(input, output, annotation_flat=None,reference_sequence=None):
        input_dep = input + ".bait_bias_summary_metrics"

        return Job(
            [input_dep],
            [output],
            [
                ['picard_convert_sequencing_artifacts_metrics', 'module_java'],
                ['picard_convert_sequencing_artifacts_metrics', 'module_picard'],
                ['picard_convert_sequencing_artifacts_metrics', 'module_R']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar ConvertSequencingArtifactToOxoG \\
 VALIDATION_STRINGENCY=SILENT  \\
 TMP_DIR={tmp_dir} \\
 INPUT_BASE={input} \\
 OUTPUT_BASE={output} \\
 REFERENCE_SEQUENCE={reference}""".format(
            tmp_dir=global_config_parser.param('picard_convert_sequencing_artifacts_metrics', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_convert_sequencing_artifacts_metrics', 'java_other_options'),
            ram=global_config_parser.param('picard_convert_sequencing_artifacts_metrics', 'ram'),
            input=input,
            output=output,
            reference=reference_sequence if reference_sequence else global_config_parser.param('picard_convert_sequencing_artifacts_metrics', 'genome_fasta'),
            )
        )

def collect_oxog_metrics(input, output, annotation_flat=None, reference_sequence=None):

        return Job(
            [input],
            [output],
            [
                ['picard_collect_sequencing_artifacts_metrics', 'module_java'],
                ['picard_collect_sequencing_artifacts_metrics', 'module_picard'],
                ['picard_collect_sequencing_artifacts_metrics', 'module_R']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar CollectOxoGMetrics \\
 VALIDATION_STRINGENCY=SILENT  \\
 TMP_DIR={tmp_dir} \\
 INPUT={input} \\
 OUTPUT={output} \\
 DB_SNP={dbsnp} \\
 REFERENCE_SEQUENCE={reference} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
            tmp_dir=global_config_parser.param('picard_collect_oxog_metrics', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_collect_oxog_metrics', 'java_other_options'),
            ram=global_config_parser.param('picard_collect_oxog_metrics', 'ram'),
            input=input,
            output=output,
            dbsnp=global_config_parser.param('picard_collect_oxog_metrics', 'known_variants'),
            reference=reference_sequence if reference_sequence else global_config_parser.param('picard_collect_oxog_metrics', 'genome_fasta'),
            max_records_in_ram=global_config_parser.param('picard_collect_oxog_metrics', 'max_records_in_ram', param_type='int')
            )
        )

def collect_gcbias_metrics(input, output, chart, summary_file, annotation_flat=None,reference_sequence=None):

        return Job(
            [input],
            [output],
            [
                ['picard_collect_gcbias_metrics', 'module_java'],
                ['picard_collect_gcbias_metrics', 'module_picard'],
                ['picard_collect_gcbias_metrics', 'module_R']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar CollectGcBiasMetrics \\
 VALIDATION_STRINGENCY=SILENT ALSO_IGNORE_DUPLICATES=TRUE \\
 TMP_DIR={tmp_dir} \\
 INPUT={input} \\
 OUTPUT={output} \\
 CHART={chart} \\
 SUMMARY_OUTPUT={summary_file} \\
 REFERENCE_SEQUENCE={reference} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
            tmp_dir=global_config_parser.param('picard_collect_gcbias_metrics', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_collect_gcbias_metrics', 'java_other_options'),
            ram=global_config_parser.param('picard_collect_gcbias_metrics', 'ram'),
            input=input,
            output=output,
            chart=chart,
            summary_file=summary_file,
            reference=reference_sequence if reference_sequence else global_config_parser.param('picard_collect_gcbias_metrics', 'genome_fasta'),
            max_records_in_ram=global_config_parser.param('picard_collect_gcbias_metrics', 'max_records_in_ram', param_type='int')
            )
        )

def fix_mate_information(input, output):

    if global_config_parser.param('fixmate', 'module_picard').split("/")[2] < "2":
        return picard.fix_mate_information(input, output)
    else:
        return Job(
            [input],
            [output],
            [
                ['fixmate', 'module_java'],
                ['fixmate', 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar FixMateInformation \\
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate \\
 TMP_DIR={tmp_dir} \\
 INPUT={input} \\
 OUTPUT={output} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
            tmp_dir=global_config_parser.param('picard_fix_mate_information', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_fix_mate_information', 'java_other_options'),
            ram=global_config_parser.param('picard_fix_mate_information', 'ram'),
            input=input,
            output=output,
            max_records_in_ram=global_config_parser.param('picard_fix_mate_information', 'max_records_in_ram', param_type='int')
            ),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
        )

def mark_duplicates(inputs, output, metrics_file, remove_duplicates="false"):

    if not isinstance(inputs, list):
        inputs=[inputs]
    if global_config_parser.param('picard_mark_duplicates', 'module_picard').split("/")[2] < "2" and global_config_parser.param('picard_mark_duplicates', 'module_gatk').split("/")[2] < "4":
        return picard.mark_duplicates(inputs, output, metrics_file, remove_duplicates)
    elif global_config_parser.param('picard_mark_duplicates', 'module_gatk').split("/")[2] > "4":
        return gatk4.mark_duplicates(inputs, output, metrics_file, remove_duplicates)
    else:
        return Job(
            inputs,
            [output, re.sub("\.([sb])am$", ".\\1ai", output), metrics_file],
            [
                ['picard_mark_duplicates', 'module_java'],
                ['picard_mark_duplicates', 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar MarkDuplicates \\
 REMOVE_DUPLICATES={remove_duplicates} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
 TMP_DIR={tmp_dir} \\
 {inputs} \\
 OUTPUT={output} \\
 METRICS_FILE={metrics_file} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram} {other_options}""".format(
            tmp_dir=global_config_parser.param('picard_mark_duplicates', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_mark_duplicates', 'java_other_options'),
            ram=global_config_parser.param('picard_mark_duplicates', 'ram'),
            remove_duplicates=remove_duplicates,
            inputs=" \\\n  ".join(["INPUT=" + str(input) for input in inputs]),
            output=output,
            metrics_file=metrics_file,
            max_records_in_ram=global_config_parser.param('picard_mark_duplicates', 'max_records_in_ram', param_type='int'),
            other_options=global_config_parser.param('picard_mark_duplicates', 'other_options', required = False) if global_config_parser.param('picard_mark_duplicates', 'other_options', required = False) else ""
            ),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
        )

def mark_duplicates_mate_cigar(inputs, output, metrics_file, remove_duplicates="false"):
    if not isinstance(inputs, list):
        inputs=[inputs]
    if global_config_parser.param('mark_duplicates_mate_cigar', 'module_gatk').split("/")[2] > "4":
        return gatk4.mark_duplicates_mate_cigar(inputs, output, metrics_file, remove_duplicates)
    else:
        return Job(
            inputs,
            [output, re.sub("\.([sb])am$", ".\\1ai", output), metrics_file],
            [
                ['mark_duplicates_mate_cigar', 'module_java'],
                ['mark_duplicates_mate_cigar', 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar \\
 MarkDuplicatesWithMateCigar \\
 REMOVE_DUPLICATES={remove_duplicates} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
 TMP_DIR={tmp_dir} \\
 {inputs} \\
 OUTPUT={output} \\
 METRICS_FILE={metrics_file} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram} {other_options}""".format(
            tmp_dir=global_config_parser.param('mark_duplicates_mate_cigar', 'tmp_dir'),
            java_other_options=global_config_parser.param('mark_duplicates_mate_cigar', 'java_other_options'),
            ram=global_config_parser.param('mark_duplicates_mate_cigar', 'ram'),
            remove_duplicates=remove_duplicates,
            inputs=" \\\n  ".join(["INPUT=" + str(input) for input in inputs]),
            output=output,
            metrics_file=metrics_file,
            max_records_in_ram=global_config_parser.param('mark_duplicates_mate_cigar', 'max_records_in_ram', param_type='int'),
            other_options=global_config_parser.param('mark_duplicates_mate_cigar', 'other_options', required=False) if global_config_parser.param('picard_mark_duplicates', 'other_options', required=False) else ""),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
        )

def mark_duplicates_mate_cigar(inputs, output, metrics_file, remove_duplicates="false"):
    if not isinstance(inputs, list):
        inputs=[inputs]
    if global_config_parser.param('picard_mark_duplicates_mate_cigar', 'module_gatk').split("/")[2] > "4":
        return gatk4.mark_duplicates_mate_cigar(inputs, output, metrics_file, remove_duplicates)
    else:
        return Job(
            inputs,
            [output, re.sub("\.([sb])am$", ".\\1ai", output), metrics_file],
            [
                ['picard_mark_duplicates_mate_cigar', 'module_java'],
                ['picard_mark_duplicates_mate_cigar', 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar \\
 MarkDuplicatesWithMateCigar \\
 REMOVE_DUPLICATES={remove_duplicates} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
 TMP_DIR={tmp_dir} \\
 {inputs} \\
 OUTPUT={output} \\
 METRICS_FILE={metrics_file} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram} {other_options}""".format(
            tmp_dir=global_config_parser.param('picard_mark_duplicates_mate_cigar', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_mark_duplicates_mate_cigar', 'java_other_options'),
            ram=global_config_parser.param('picard_mark_duplicates_mate_cigar', 'ram'),
            remove_duplicates=remove_duplicates,
            inputs=" \\\n  ".join(["INPUT=" + str(input) for input in inputs]),
            output=output,
            metrics_file=metrics_file,
            max_records_in_ram=global_config_parser.param('picard_mark_duplicates_mate_cigar', 'max_records_in_ram', param_type='int'),
            other_options=global_config_parser.param('picard_mark_duplicates_mate_cigar', 'other_options', required = False) if global_config_parser.param('picard_mark_duplicates', 'other_options', required = False) else ""),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
        )

def merge_sam_files(inputs, output):

    if not isinstance(inputs, list):
        inputs=[inputs]
    if global_config_parser.param('picard_merge_sam_files', 'module_picard').split("/")[2] < "2":
        return picard.merge_sam_files(inputs, output)
    else:
        return Job(
            inputs,
            [output, re.sub("\.([sb])am$", ".\\1ai", output)],
            [
                ['picard_merge_sam_files', 'module_java'],
                ['picard_merge_sam_files', 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar MergeSamFiles \\
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \\
 TMP_DIR={tmp_dir} \\
 {inputs} \\
 OUTPUT={output} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
            tmp_dir=global_config_parser.param('picard_merge_sam_files', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_merge_sam_files', 'java_other_options'),
            ram=global_config_parser.param('picard_merge_sam_files', 'ram'),
            inputs=" \\\n ".join(["INPUT=" + input for input in inputs]),
            output=output,
            max_records_in_ram=global_config_parser.param('picard_merge_sam_files', 'max_records_in_ram', param_type='int')
            ),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output)]
        )

# Reorder BAM/SAM files based on reference/dictionary
def reorder_sam(input, output):

    
    if global_config_parser.param('reorder_sam', 'module_picard').split("/")[2] < "2":
        return picard.reorder_sam(input, output)
    else:
        return Job(
            [input],
            [output],
            [
                ['reorder_sam', 'module_java'],
                ['reorder_sam', 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar ReorderSam \\
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
 TMP_DIR={tmp_dir} \\
 INPUT={input} \\
 OUTPUT={output} \\
 REFERENCE={reference} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
            tmp_dir=global_config_parser.param('picard_reorder_sam', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_reorder_sam', 'java_other_options'),
            ram=global_config_parser.param('picard_reorder_sam', 'ram'),
            input=input,
            output=output,
            reference=global_config_parser.param('picard_reorder_sam', 'genome_fasta', param_type='filepath'),
            max_records_in_ram=global_config_parser.param('picard_reorder_sam', 'max_records_in_ram', param_type='int')
            ),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output)]
        )

# Convert SAM/BAM file to fastq format
def sam_to_fastq(input, fastq, second_end_fastq=None):

    if global_config_parser.param('picard_sam_to_fastq', 'module_picard').split("/")[2] < "2":
        return picard.sam_to_fastq(input, fastq, second_end_fastq)
    else:
        return Job(
            [input],
            [fastq, second_end_fastq],
            [
                ['picard_sam_to_fastq', 'module_java'],
                ['picard_sam_to_fastq', 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar SamToFastq \\
 VALIDATION_STRINGENCY=LENIENT \\
 CREATE_MD5_FILE=TRUE \\
 INPUT={input} \\
 FASTQ={fastq}{second_end_fastq}""".format(
            tmp_dir=global_config_parser.param('picard_sam_to_fastq', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_sam_to_fastq', 'java_other_options'),
            ram=global_config_parser.param('picard_sam_to_fastq', 'ram'),
            input=input,
            fastq=fastq,
            second_end_fastq=" \\\n  SECOND_END_FASTQ=" + second_end_fastq if second_end_fastq else ""
            ),
            removable_files=[fastq, second_end_fastq]
        )

def sort_sam(input, output, sort_order="coordinate", ini_section='picard_sort_sam'):

    if global_config_parser.param(ini_section, 'module_picard').split("/")[2] < "2":
        return picard.sort_sam(input, output, sort_order, ini_section)
    else:
        return Job(
            [input],
            # Add SAM/BAM index as output only when writing a coordinate-sorted BAM file
            [output, re.sub("\.([sb])am$", ".\\1ai", output) if sort_order == "coordinate" else None],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_picard']
            ],
            command="""\
 java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar SortSam \\
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
 TMP_DIR={tmp_dir} \\
 INPUT={input} \\
 OUTPUT={output} \\
 SORT_ORDER={sort_order} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
            tmp_dir=global_config_parser.param(ini_section, 'tmp_dir'),
            java_other_options=global_config_parser.param(ini_section, 'java_other_options'),
            ram=global_config_parser.param(ini_section, 'ram'),
            input=input,
            output=output,
            sort_order=sort_order,
            max_records_in_ram=global_config_parser.param(ini_section, 'max_records_in_ram', param_type='int')
            ),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output) if sort_order == "coordinate" else None]
        )

def sort_vcfs(inputs, output, ini_section='picard_sort_vcf'):

    if not isinstance(inputs, list):
        inputs=[inputs]
    if global_config_parser.param(ini_section, 'module_picard').split("/")[2] < "2":
        return picard.sort_vcfs(inputs, output, ini_section)
    else:
        return Job(
            inputs,
            # Add SAM/BAM index as output only when writing a coordinate-sorted BAM file
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar SortVcf \\
 VALIDATION_STRINGENCY=SILENT \\
 TMP_DIR={tmp_dir} \\
 {inputs} \\
 OUTPUT={output} \\
 SEQUENCE_DICTIONARY={seq_dict}""".format(
            tmp_dir=global_config_parser.param(ini_section, 'tmp_dir'),
            java_other_options=global_config_parser.param(ini_section, 'java_other_options'),
            ram=global_config_parser.param(ini_section, 'ram'),
            inputs=" \\\n  ".join(["INPUT=" + input for input in inputs]),
            output=output,
            seq_dict=global_config_parser.param(ini_section, 'genome_dictionary', param_type='filepath')
            )
        )

def mergeVcfs(variants, output):
    return Job(
            variants,
            [output],
            [
                ['picard_merge_vcfs', 'module_java'],
                ['picard_merge_vcfs', 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar \\
  MergeVcfs {options} \\
  --TMP_DIR {tmp_dir} \\
  --REFERENCE_SEQUENCE {reference}{variants} \\
  --OUTPUT {output} \\
  --MAX_RECORDS_IN_RAM {max_records_in_ram}""".format(
                tmp_dir=global_config_parser.param('picard_merge_vcfs', 'tmp_dir'),
                java_other_options=global_config_parser.param('picard_merge_vcfs', 'java_other_options'),
                ram=global_config_parser.param('picard_merge_vcfs', 'ram'),
                options=global_config_parser.param('picard_merge_vcfs', 'options'),
                reference=global_config_parser.param('picard_merge_vcfs', 'genome_fasta', param_type='filepath'),
                variants="".join(" \\\n  --INPUT " + variant for variant in variants),
                output=output,
                max_records_in_ram=global_config_parser.param('picard_merge_vcfs', 'max_records_in_ram', param_type='int')
        )
  )

def collect_rna_metrics(input, output, annotation_flat=None,reference_sequence=None):

    if global_config_parser.param('picard_collect_rna_metrics', 'module_picard').split("/")[2] < "2":
        return picard.collect_rna_metrics(input, output, annotation_flat, reference_sequence)
    else:
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
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \\
 VALIDATION_STRINGENCY=SILENT  \\
 TMP_DIR={tmp_dir} \\
 INPUT={input} \\
 OUTPUT={output} \\
 REF_FLAT={ref_flat} \\
 STRAND_SPECIFICITY={strand_specificity} \\
 MINIMUM_LENGTH={min_length} \\
 REFERENCE_SEQUENCE={reference} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
            tmp_dir=global_config_parser.param('picard_collect_rna_metrics', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_collect_rna_metrics', 'java_other_options'),
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

def add_read_groups(input, output, readgroup, library, processing_unit, sample, sort_order="coordinate"):
    if global_config_parser.param('picard_add_read_groups', 'module_picard').split("/")[2] < "2":
        return picard.add_read_groups(input, output, readgroup, library, processing_unit, sample, sort_order)
    else:
        return Job(
            [input],
            # collect specific RNA metrics (exon rate, strand specificity, etc...)
            [output, re.sub("\.([sb])am$", ".\\1ai", output)],
            [
                ['picard_add_read_groups', 'module_java'],
                ['picard_add_read_groups', 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups \\
 VALIDATION_STRINGENCY=SILENT  \\
 TMP_DIR={tmp_dir} \\
 CREATE_INDEX=true \\
 INPUT={input} \\
 OUTPUT={output} \\
 SORT_ORDER=\"{sort_order}\" \\
 RGID=\"{readgroup}\" \\
 RGLB=\"{library}\" \\
 RGPL=\"{platform}\" \\
 RGPU=\"run{processing_unit}\" \\
 RGSM=\"{sample}\" \\
 {sequencing_center}""".format(
                tmp_dir=global_config_parser.param('picard_add_read_groups', 'tmp_dir'),
                java_other_options=global_config_parser.param('picard_add_read_groups', 'java_other_options'),
                ram=global_config_parser.param('picard_add_read_groups', 'ram'),
                input=input,
                output=output,
                sort_order=sort_order,
                readgroup=readgroup,
                library=library,
                platform=global_config_parser.param('picard_add_read_groups', 'platform'),
                processing_unit=processing_unit,
                sample=sample,
                sequencing_center=("RGCN=\"" + global_config_parser.param(
                    'picard_add_read_groups', 'sequencing_center') + "\""
                                   if global_config_parser.param(
                    'picard_add_read_groups', 'sequencing_center', required=False) else "")
            )
        )

def bed2interval_list(
    dictionary,
    bed,
    output
    ):
    
    if global_config_parser.param('picard_bed2interval_list', 'module_picard').split("/")[2] < "2":
        return picard.bed2interval_list(
            dictionary,
            bed,
            output
        )
    return Job(
        [dictionary, bed],
        [output],
        [
            ['picard_bed2interval_list', 'module_java'],
            ['picard_bed2interval_list', 'module_picard']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar BedToIntervalList \\
  INPUT={bed} \\
  SEQUENCE_DICTIONARY={dictionary} \\
  OUTPUT={output}""".format(
            tmp_dir=global_config_parser.param('picard_bed2interval_list', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_bed2interval_list', 'java_other_options'),
            ram=global_config_parser.param('picard_bed2interval_list', 'ram'),
            dictionary=dictionary if dictionary else global_config_parser.param('picard_bed2interval_list', 'genome_dictionary', param_type='filepath'),
            bed=bed,
            output=output,
        )
    )

def interval_list2bed(input, output):
    return Job(
        [input],
        [output],
        [
            ['picard_interval_list2bed', 'module_java'],
            ['picard_interval_list2bed', 'module_picard']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar IntervalListToBed \\
  INPUT={input} \\
  OUTPUT={output}""".format(
            tmp_dir=global_config_parser.param('picard_interval_list2bed', 'tmp_dir'),
            java_other_options=global_config_parser.param('picard_interval_list2bed', 'java_other_options'),
            ram=global_config_parser.param('picard_interval_list2bed', 'ram'),
            input=input,
            output=output
            )
        )

def scatterIntervalsByNs(reference,
                  output,
                  options= None):
#                  exclude_intervals=None):

    return Job(
        [reference],
        [output],
        [
            ['picard_ScatterIntervalsByNs', 'module_java'],
            ['picard_ScatterIntervalsByNs', 'module_picard']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/picard.jar \\
  ScatterIntervalsByNs {options} \\
  REFERENCE={reference} \\
  OUTPUT={output}""".format(
            tmp_dir=global_config_parser.param('picard_ScatterIntervalsByNs', 'tmp_dir'),
            options=options if options else global_config_parser.param('picard_ScatterIntervalsByNs', 'options'),
            java_other_options=global_config_parser.param('picard_ScatterIntervalsByNs', 'java_other_options'),
            ram=global_config_parser.param('picard_ScatterIntervalsByNs', 'ram'),
            reference=reference if reference else global_config_parser.param('picard_ScatterIntervalsByNs', 'genome_fasta', param_type='filepath'),
#            exclude_intervals=exclude_intervals if exclude_intervals else "".join(" \\\n  --excludeIntervals " + exclude_interval for exclude_interval in exclude_intervals),
            output=output
        )
    )
