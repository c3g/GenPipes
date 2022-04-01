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

# MUGQIC Modules
from core.config import config
from core.job import Job
from . import picard2

def build_bam_index(input, output):

    if config.param('build_bam_index', 'module_picard').split("/")[2] >= "2":
        return picard2.build_bam_index(input, output)
    else:
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

    if config.param('picard_calculate_hs_metrics', 'module_picard').split("/")[2] >= "2":
        return picard2.calculate_hs_metrics(input, output, intervals, reference_sequence)
    else:
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
            reference_sequence=reference_sequence if reference_sequence else config.param('picard_calculate_hs_metrics', 'genome_fasta', param_type='filepath')
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

    if config.param('picard_collect_multiple_metrics', 'module_picard').split("/")[2] >= "2":
        return picard2.collect_multiple_metrics(input, output, reference_sequence, library_type)
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
            reference_sequence=reference_sequence if reference_sequence else config.param('picard_collect_multiple_metrics', 'genome_fasta', param_type='filepath'),
            input=input,
            output=output,
            max_records_in_ram=config.param('picard_collect_multiple_metrics', 'max_records_in_ram', param_type='int')
            )
        )

def fix_mate_information(input, output):

    if config.param('fixmate', 'module_picard').split("/")[2] >= "2":
        return picard2.fix_mate_information(input, output)
    else:
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
            max_records_in_ram=config.param('picard_fix_mate_information', 'max_records_in_ram', param_type='int')
            ),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
        )

def mark_duplicates(inputs, output, metrics_file, remove_duplicates="false"):

    if not isinstance(inputs, list):
        inputs=[inputs]
    if config.param('picard_mark_duplicates', 'module_picard').split("/")[2] >= "2":
        return picard2.mark_duplicates(inputs, output, metrics_file, remove_duplicates)
    else:
        return Job(
            inputs,
            [output, re.sub("\.([sb])am$", ".\\1ai", output), metrics_file],
            [
                ['picard_mark_duplicates', 'module_java'],
                ['picard_mark_duplicates', 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/MarkDuplicates.jar \\
 REMOVE_DUPLICATES={remove_duplicates} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
 TMP_DIR={tmp_dir} \\
 {inputs} \\
 OUTPUT={output} \\
 METRICS_FILE={metrics_file} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram} {other_options}""".format(
            tmp_dir=config.param('picard_mark_duplicates', 'tmp_dir'),
            java_other_options=config.param('picard_mark_duplicates', 'java_other_options'),
            ram=config.param('picard_mark_duplicates', 'ram'),
            remove_duplicates=remove_duplicates,
            inputs=" \\\n  ".join(["INPUT=" + input for input in inputs]),
            output=output,
            metrics_file=metrics_file,
            max_records_in_ram=config.param('picard_mark_duplicates', 'max_records_in_ram', param_type='int'),
            other_options= config.param('picard_mark_duplicates', 'other_options',required = False) if config.param('picard_mark_duplicates', 'other_options',required = False) else ""
            ),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output), output + ".md5"]
        )

def merge_sam_files(inputs, output):

    if not isinstance(inputs, list):
        inputs=[inputs]
    if config.param('picard_merge_sam_files', 'module_picard').split("/")[2] >= "2":
        return picard2.merge_sam_files(inputs, output)
    else:
        return Job(
            inputs,
            [output, re.sub("\.([sb])am$", ".\\1ai", output)],
            [
                ['picard_merge_sam_files', 'module_java'],
                ['picard_merge_sam_files', 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/MergeSamFiles.jar \\
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true USE_THREADING=true \\
  TMP_DIR={tmp_dir} \\
  {inputs} \\
  OUTPUT={output} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=config.param('picard_merge_sam_files', 'tmp_dir'),
        java_other_options=config.param('picard_merge_sam_files', 'java_other_options'),
        ram=config.param('picard_merge_sam_files', 'ram'),
        inputs=" \\\n  ".join(["INPUT=" + input for input in inputs]),
        output=output,
        max_records_in_ram=config.param('picard_merge_sam_files', 'max_records_in_ram', param_type='int')
        ),
        removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output)]
    )

# Reorder BAM/SAM files based on reference/dictionary
def reorder_sam(input, output):

    if config.param('reorder_sam', 'module_picard').split("/")[2] >= "2":
        return picard2.reorder_sam(input, output)
    else:
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
            reference=config.param('picard_reorder_sam', 'genome_fasta', param_type='filepath'),
            max_records_in_ram=config.param('picard_reorder_sam', 'max_records_in_ram', param_type='int')
            ),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output)]
        )

# Convert SAM/BAM file to fastq format
def sam_to_fastq(input, fastq, second_end_fastq=None):

    if config.param('picard_sam_to_fastq', 'module_picard').split("/")[2] >= "2":
        return picard2.sam_to_fastq(input, fastq, second_end_fastq)
    else:
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
 CREATE_MD5_FILE=TRUE \\
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

def sort_sam(input, output, sort_order="coordinate", ini_section='picard_sort_sam'):

    if config.param(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.sort_sam(input, output, sort_order, ini_section)
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
            max_records_in_ram=config.param(ini_section, 'max_records_in_ram', param_type='int')
            ),
            removable_files=[output, re.sub("\.([sb])am$", ".\\1ai", output) if sort_order == "coordinate" else None]
        )

def sort_vcfs(inputs, output, ini_section='picard_sort_vcf'):

    if not isinstance(inputs, list):
        inputs=[inputs]
    if config.param(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.sort_vcfs(inputs, output, ini_section)
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
            seq_dict=config.param(ini_section, 'genome_dictionary', param_type='filepath')
            )
        )

def collect_rna_metrics(input, output, annotation_flat=None,reference_sequence=None):

    if config.param('picard_collect_rna_metrics', 'module_picard').split("/")[2] >= "2":
        return picard2.collect_rna_metrics(input, output, annotation_flat,reference_sequence)
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
            min_length=config.param('picard_collect_rna_metrics', 'minimum_length', param_type='int'),
            reference=reference_sequence if reference_sequence else config.param('picard_collect_rna_metrics', 'genome_fasta'),
            max_records_in_ram=config.param('picard_collect_rna_metrics', 'max_records_in_ram', param_type='int')
            )
        )


def add_read_groups(input, output, readgroup, library, processing_unit, sample, sort_order="coordinate"):

    if config.param('picard_add_read_groups', 'module_picard').split("/")[2] >= "2":
        return picard2.add_read_groups(input, output, readgroup, library, processing_unit, sample, sort_order)
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
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/AddOrReplaceReadGroups.jar \\
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
                tmp_dir=config.param('picard_add_read_groups', 'tmp_dir'),
                java_other_options=config.param('picard_add_read_groups', 'java_other_options'),
                ram=config.param('picard_add_read_groups', 'ram'),
                input=input,
                output=output,
                sort_order=sort_order,
                readgroup=readgroup,
                library=library,
                platform=config.param('picard_add_read_groups', 'platform'),
                processing_unit=processing_unit,
                sample=sample,
                sequencing_center=("RGCN=\"" + config.param('picard_add_read_groups', 'sequencing_center') + "\"") if config.param('picard_add_read_groups', 'sequencing_center', required=False) else ""
            )
        )


def bed2interval_list(
    dictionary,
    bed,
    output
    ):

    if config.param('picard_bed2interval_list', 'module_picard').split("/")[2] >= "2":
        return picard2.bed2interval_list(
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
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/BedToIntervalList.jar \\
  INPUT={bed} \\
  SEQUENCE_DICTIONARY={dictionary} \\
  OUTPUT={output}""".format(
            tmp_dir=config.param('picard_bed2interval_list', 'tmp_dir'),
            java_other_options=config.param('picard_bed2interval_list', 'java_other_options'),
            ram=config.param('picard_bed2interval_list', 'ram'),
            dictionary=dictionary if dictionary else config.param('picard_bed2interval_list', 'genome_dictionary', param_type='filepath'),
            bed=bed,
            output=output
            )
        )
