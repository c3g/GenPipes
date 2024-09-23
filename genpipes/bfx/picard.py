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

# GenPipes Modules
from ..core.config import global_conf
from ..core.job import Job
from . import picard2

def build_bam_index(input, output, ini_section='build_bam_index'):

    if global_conf.global_get(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.build_bam_index(input, output)
    else:
        return Job(
            [input],
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/BuildBamIndex.jar \\
 VALIDATION_STRINGENCY=SILENT \\
 INPUT={input} \\
 OUTPUT={output} """.format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            input=input,
            output=output,
            )
        )

def calculate_hs_metrics(input, output, intervals, reference_sequence=None, ini_section='picard_calculate_hs_metrics'):

    baits_intervals = ""
    baits_intervals = global_conf.global_get(ini_section, 'baits_intervals', required = False)

    if global_conf.global_get(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.calculate_hs_metrics(input, output, intervals, reference_sequence, ini_section=ini_section)
    else:
        return Job(
            [input, intervals],
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/CalculateHsMetrics.jar \\
 TMP_DIR={tmp_dir} \\
 INPUT={input} \\
 OUTPUT={output} \\
 BAIT_INTERVALS={baits} \\
 TARGET_INTERVALS={intervals} \\
 REFERENCE_SEQUENCE={reference_sequence}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            input=input,
            output=output,
            intervals=intervals,
            baits=baits_intervals if baits_intervals != "" else intervals,
            reference_sequence=reference_sequence if reference_sequence else global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath')
            )
        )

def collect_multiple_metrics(input, output, reference_sequence=None , library_type="PAIRED_END", ini_section='picard_collect_multiple_metrics'):

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

    if global_conf.global_get(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.collect_multiple_metrics(input, output, reference_sequence, library_type, ini_section=ini_section)
    else:
        return Job(
            [input],
            outputs,
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_picard'],
                [ini_section, 'module_R']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/CollectMultipleMetrics.jar \\
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \\
 TMP_DIR={tmp_dir} \\
 REFERENCE_SEQUENCE={reference_sequence} \\
 INPUT={input} \\
 OUTPUT={output} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                reference_sequence=reference_sequence if reference_sequence else global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
                input=input,
                output=output,
                max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            ),
            report_files=outputs
        )

def fix_mate_information(input, output, ini_section='picard_fix_mate_information'):

    if global_conf.global_get(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.fix_mate_information(input, output, ini_section=ini_section)
    else:
        return Job(
            [input],
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/FixMateInformation.jar \\
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate \\
 TMP_DIR={tmp_dir} \\
 INPUT={input} \\
 OUTPUT={output} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            input=input,
            output=output,
            max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            ),
            removable_files=[output, re.sub(r"\.([sb])am$", ".\\1ai", output), output + ".md5"]
        )

def mark_duplicates(
    inputs,
    output,
    metrics_file,
    remove_duplicates="false",
    create_index=True,
    ini_section='picard_mark_duplicates'
    ):

    if not isinstance(inputs, list):
        inputs=[inputs]
    if global_conf.global_get(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.mark_duplicates(inputs, output, metrics_file, remove_duplicates, ini_section=ini_section)
    else:
        outputs = [output, metrics_file]
        if create_index:
            outputs.append(re.sub(r"\.([sb])am$", ".\\1ai", output))
        return Job(
            inputs,
            outputs,
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_picard']
            ],
            command="""\
rm -rf {output}.part && \\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/MarkDuplicates.jar \\
 REMOVE_DUPLICATES={remove_duplicates} VALIDATION_STRINGENCY=SILENT {create_index} \\
 TMP_DIR={tmp_dir} \\
 {inputs} \\
 OUTPUT={output} \\
 METRICS_FILE={metrics_file} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram} {other_options}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                remove_duplicates=remove_duplicates,
                create_index="CREATE_INDEX=true" if create_index else "",
                inputs=" \\\n  ".join(["INPUT=" + input for input in inputs]),
                output=output,
                metrics_file=metrics_file,
                max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int'),
                other_options= global_conf.global_get(ini_section, 'other_options',required = False)
            ),
            removable_files=[output, re.sub(r"\.([sb])am$", ".\\1ai", output), output + ".md5"]
        )

def merge_sam_files(inputs, output, ini_section='picard_merge_sam_files'):

    if not isinstance(inputs, list):
        inputs=[inputs]
    if global_conf.global_get(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.merge_sam_files(inputs, output, ini_section=ini_section)
    else:
        return Job(
            inputs,
            [output, re.sub(r"\.([sb])am$", ".\\1ai", output)],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/MergeSamFiles.jar \\
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true USE_THREADING=true \\
  TMP_DIR={tmp_dir} \\
  {inputs} \\
  OUTPUT={output} \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
        tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
        java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
        ram=global_conf.global_get(ini_section, 'ram'),
        inputs=" \\\n  ".join(["INPUT=" + input for input in inputs]),
        output=output,
        max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
        ),
        removable_files=[output, re.sub(r"\.([sb])am$", ".\\1ai", output)]
    )

# Reorder BAM/SAM files based on reference/dictionary
def reorder_sam(input, output, ini_section='picard_reorder_sam'):

    if global_conf.global_get(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.reorder_sam(input, output, ini_section=ini_section)
    else:
        return Job(
            [input],
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/ReorderSam.jar \\
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
 TMP_DIR={tmp_dir} \\
 INPUT={input} \\
 OUTPUT={output} \\
 REFERENCE={reference} \\
 MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            input=input,
            output=output,
            reference=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
            max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            ),
            removable_files=[output, re.sub(r"\.([sb])am$", ".\\1ai", output)]
        )

# Convert SAM/BAM file to fastq format
def sam_to_fastq(input, fastq, second_end_fastq=None, ini_section='picard_sam_to_fastq'):

    if global_conf.global_get(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.sam_to_fastq(input, fastq, second_end_fastq, ini_section=ini_section)
    else:
        return Job(
            [input],
            [fastq, second_end_fastq],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_picard']
            ],
            command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/SamToFastq.jar \\
  {other_options} \\
  CREATE_MD5_FILE=TRUE \\
  INPUT={input} \\
  FASTQ={fastq}{second_end_fastq}""".format(
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                other_options=global_conf.global_get(ini_section, 'other_options'),
                input=input,
                fastq=fastq,
                second_end_fastq=" \\\n  SECOND_END_FASTQ=" + second_end_fastq if second_end_fastq else ""
            ),
            removable_files=[fastq, second_end_fastq]
        )

def sort_sam(input, output, sort_order="coordinate", ini_section='picard_sort_sam'):

    if global_conf.global_get(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.sort_sam(input, output, sort_order, ini_section)
    else:
        return Job(
            [input],
            # Add SAM/BAM index as output only when writing a coordinate-sorted BAM file
            [output, re.sub(r"\.([sb])am$", ".\\1ai", output) if sort_order == "coordinate" else None],
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
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            input=input,
            output=output,
            sort_order=sort_order,
            max_records_in_ram=global_conf.global_get(ini_section, 'max_records_in_ram', param_type='int')
            ),
            removable_files=[output, re.sub(r"\.([sb])am$", ".\\1ai", output) if sort_order == "coordinate" else None]
        )

def sort_vcfs(inputs, output, ini_section='picard_sort_vcf'):

    if not isinstance(inputs, list):
        inputs=[inputs]
    if global_conf.global_get(ini_section, 'module_picard').split("/")[2] >= "2":
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
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            inputs=" \\\n  ".join(["INPUT=" + input for input in inputs]),
            output=output,
            seq_dict=global_conf.global_get(ini_section, 'genome_dictionary', param_type='filepath')
            )
        )

def collect_rna_metrics(input, output, annotation_flat=None,reference_sequence=None, ini_section='picard_collect_rna_metrics'):

    if global_conf.global_get(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.collect_rna_metrics(input, output, annotation_flat,reference_sequence, ini_section=ini_section)
    else:
        return Job(
            [input],
            # collect specific RNA metrics (exon rate, strand specificity, etc...)
            [output],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_picard'],
                [ini_section, 'module_R']
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
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
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


def add_read_groups(input, output, readgroup, library, processing_unit, sample, sort_order="coordinate", ini_section='picard_add_read_groups'):

    if global_conf.global_get(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.add_read_groups(input, output, readgroup, library, processing_unit, sample, sort_order, ini_section=ini_section)
    else:
        return Job(
            [input],
            # collect specific RNA metrics (exon rate, strand specificity, etc...)
            [output, re.sub(r"\.([sb])am$", ".\\1ai", output)],
            [
                [ini_section, 'module_java'],
                [ini_section, 'module_picard']
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
                tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
                java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
                ram=global_conf.global_get(ini_section, 'ram'),
                input=input,
                output=output,
                sort_order=sort_order,
                readgroup=readgroup,
                library=library,
                platform=global_conf.global_get(ini_section, 'platform'),
                processing_unit=processing_unit,
                sample=sample,
                sequencing_center=("RGCN=\"" + global_conf.global_get(ini_section, 'sequencing_center') + "\"") if global_conf.global_get(ini_section, 'sequencing_center', required=False) else ""
            )
        )


def bed2interval_list(
    dictionary,
    bed,
    output,
    ini_section='picard_bed2interval_list'
    ):

    if global_conf.global_get(ini_section, 'module_picard').split("/")[2] >= "2":
        return picard2.bed2interval_list(
            dictionary,
            bed,
            output,
            ini_section=ini_section
        )

    return Job(
        [dictionary, bed],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_picard']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PICARD_HOME/BedToIntervalList.jar \\
  INPUT={bed} \\
  SEQUENCE_DICTIONARY={dictionary} \\
  OUTPUT={output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            dictionary=dictionary if dictionary else global_conf.global_get(ini_section, 'genome_dictionary', param_type='filepath'),
            bed=bed,
            output=output
            )
        )
