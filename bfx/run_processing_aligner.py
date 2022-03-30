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
import os
import re
import sys

# MUGQIC Modules
from core.job import Job, concat_jobs, pipe_jobs
from core.config import config
from bfx import bvatools
from bfx import verify_bam_id
from bfx import bwa
from bfx import metrics
from bfx import picard
from bfx import star
from bfx import cellranger
from bfx import tools
from bfx import bash_cmd as bash

import logging
log = logging.getLogger(__name__)

class RunProcessingAligner(object):
    def __init__(self, output_dir, genome_folder, platform):
        self._output_dir = output_dir
        self._genome_folder = genome_folder
        self._platform = platform

    @property
    def output_dir(self):
        return self._output_dir

    @property
    def genome_folder(self):
        return self._genome_folder

    @property
    def platform(self):
        return self._platform

    def get_reference_index(self):
        raise NotImplementedError("Please Implement this method")

    def get_reference_file(self):
        raise NotImplementedError("Please Implement this method")

    def get_alignment_job(self, readset):
        raise NotImplementedError("Please Implement this method")

    def get_metrics_jobs(self, readset):
        raise NotImplementedError("Please Implement this method")

    def get_annotation_files(self):
        raise NotImplementedError("Please Implement this method")

    @staticmethod
    def get_rg_tag(readset, platform, ini_section):
        return "'@RG" + \
               "\tID:" + readset.library + "_" + readset.run + "_" + readset.lane + \
               "\tSM:" + readset.sample.name + \
               "\tLB:" + readset.library + \
               "\tPU:run" + readset.run + "_" + readset.lane + \
               ("\tCN:" + config.param(ini_section, 'sequencing_center')
                if config.param(ini_section, 'sequencing_center', required=False) else "") + \
               "\tPL:" + platform + \
               "'"

    def verify_bam_id(self, readset, type="DNA"):
        """
            verifyBamID is a software that verifies whether the reads in particular file match previously known
            genotypes for an individual (or group of individuals), and checks whether the reads are contaminated
            as a mixture of two samples. verifyBamID can detect sample contamination and swaps when external
            genotypes are available. When external genotypes are not available, verifyBamID still robustly
            detects sample swaps.
        """
        jobs = []
        known_variants_annotated_filtered = ""
        if type == "RNA":
            if len(readset.annotation_files) > 3 and os.path.isfile(readset.annotation_files[3]):
                known_variants_annotated = readset.annotation_files[3]
            else:
                return jobs
        else:
            if len(readset.annotation_files) > 0 and os.path.isfile(readset.annotation_files[0]):
                known_variants_annotated = readset.annotation_files[0]
            else:
                return jobs

        known_variants_annotated_filtered = known_variants_annotated

        if known_variants_annotated_filtered:

            input_bam = readset.bam + ".bam"
            output_prefix = readset.bam + ".metrics.verifyBamId"
            jobs.append(
                concat_jobs([
                    verify_bam_id.verify(
                        input_bam,
                        output_prefix,
                        vcf=known_variants_annotated_filtered,
                    ),
                    # the first column starts with a # (a comment for nanuq) so we remove the column and output the result
                    # in a file with the name supported by nanuq
                    Job(
                        [output_prefix + ".selfSM"],
                        [output_prefix + ".tsv"],
                        command="cut -f2- " + output_prefix + ".selfSM > " + output_prefix + ".tsv"
                )],
                name="verify_bam_id." + readset.name + "." + readset.run + "." + readset.lane,
                report_files=[readset.bam + ".metrics.verifyBamId.selfSM"],
                samples=[readset.sample]
            ))

        return jobs

    @staticmethod
    def find_10x_synonym_reference(parent_genome_folder, synonym, option):
        trans_path = ""
        assembly_name = os.path.basename(parent_genome_folder)
        species = assembly_name.split(".")[0]
        main_genome_folder = os.path.dirname(parent_genome_folder)
        new_assembly = species + "." + synonym
        synonym_path = os.path.join(main_genome_folder, new_assembly)
        ini_file = os.path.join(synonym_path + os.sep + new_assembly + ".ini")

        if os.path.isfile(ini_file):
            config.parse_files([ini_file])

            if config.has_option("DEFAULT", option):
                transcriptome = config.get("DEFAULT", option)
                if transcriptome:
                    trans_path = os.path.join(synonym_path, "genome", "10xGenomics", transcriptome)
        return trans_path

class NullRunProcessingAligner(RunProcessingAligner):
    """ A no-op aligner used when we want to skip the alignment """

    def get_reference_index(self):
        return None

    def get_reference_file(self):
        return None

    def get_alignment_job(self, readset):
        return []

    def get_metrics_jobs(self, readset):
        return []

    def get_annotation_files(self):
        return []

class BwaRunProcessingAligner(RunProcessingAligner):
    downloaded_bed_files = []
    created_interval_lists = []

    def get_reference_index(self):
        folder_name = os.path.basename(self.genome_folder)
        return os.path.join(
            self.genome_folder,
            "genome",
            "bwa_index",
            folder_name + ".fa"
        )

    def get_reference_file(self):
        folder_name = os.path.basename(self.genome_folder)
        return os.path.join(
            self.genome_folder,
            "genome",
            folder_name + ".fa"
        )

    def get_annotation_files(self):
        folder_name = os.path.basename(self.genome_folder)
        ini_file = os.path.join(self.genome_folder + os.sep + folder_name + ".ini")
        if os.path.isfile(ini_file):
            config.parse_files([ini_file])

            section = "DEFAULT"
            dbsnp_option_name = "dbsnp_version"
            af_option_name = "population_AF"

            if config.has_option(section, dbsnp_option_name) and config.has_option(section, af_option_name):
                dbsnp_version = config.get(section, dbsnp_option_name)
                af_name = config.get(section, af_option_name)
                return [
                    os.path.join(
                        self.genome_folder,
                        "annotations",
                        folder_name + ".dbSNP" + dbsnp_version + "_" + af_name + ".vcf.gz"
                    )
                ]

        return []

    def get_alignment_job(self, readset):
        output = readset.bam + ".bam"
        job = concat_jobs([
            bash.mkdir(os.path.dirname(output)),
            pipe_jobs([
                bwa.mem(
                    readset.fastq1,
                    readset.fastq2,
                    read_group=RunProcessingAligner.get_rg_tag(readset, self.platform, 'bwa_mem'),
                    ref=readset.aligner_reference_index
                ),
                picard.sort_sam(
                    "/dev/stdin",
                    output,
                    "coordinate"
                )]
            )],
            name="bwa_mem_picard_sort_sam." + readset.name + "." + readset.run + "." + readset.lane,
            samples=[readset.sample]
        )

        return job

    def get_metrics_jobs(self, readset):
        jobs = []

        input_file_prefix = readset.bam + '.'
        input = input_file_prefix + "bam"

        job = picard.collect_multiple_metrics(
            input,
            input_file_prefix + "metrics",
            reference_sequence=readset.reference_file
        )
        job.name = "picard_collect_multiple_metrics." + readset.name + ".met" + "." + readset.run + "." + readset.lane
        job.samples = [readset.sample]
        job.report_files = [
            input_file_prefix + "metrics.alignment_summary_metrics",
            input_file_prefix + "metrics.insert_size_metrics"
        ]
        jobs.append(job)

        if readset.beds:
            coverage_bed = readset.beds[0]
            full_coverage_bed = os.path.join(config.param('DEFAULT', 'bed_path', param_type='dirpath', required=False), coverage_bed)
            if not os.path.isfile(full_coverage_bed):
                full_coverage_bed = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", "bed", coverage_bed)
        else:
            coverage_bed = None
            full_coverage_bed = None

        if coverage_bed:
            if (not os.path.exists(full_coverage_bed)) and (coverage_bed not in BwaRunProcessingAligner.downloaded_bed_files):
                # Download the bed file
                command = config.param('DEFAULT', 'fetch_bed_file_command').format(
                    output_directory=self.output_dir,
                    filename=coverage_bed
                )
                job = Job([], [full_coverage_bed], command=command, name="bed_download." + coverage_bed, samples=[readset.sample])
                BwaRunProcessingAligner.downloaded_bed_files.append(coverage_bed)
                jobs.append(job)

            interval_list = re.sub("\.[^.]+$", ".interval_list", coverage_bed)

            if interval_list not in BwaRunProcessingAligner.created_interval_lists:
                # Create one job to generate the interval list from the bed file
                ref_dict = os.path.splitext(readset.reference_file)[0] + '.dict'
                job = tools.bed2interval_list(
                    full_coverage_bed,
                    interval_list,
                    ref_dict
                )
                job.name = "interval_list." + coverage_bed
                job.samples = [readset.sample]
                BwaRunProcessingAligner.created_interval_lists.append(interval_list)
                jobs.append(job)

            job = picard.calculate_hs_metrics(input_file_prefix + "bam", input_file_prefix + "metrics.onTarget.txt",
                                              interval_list, reference_sequence=readset.reference_file)
            job.name = "picard_calculate_hs_metrics." + readset.name + ".hs" + "." + readset.run + "." + readset.lane
            job.samples = [readset.sample]
            jobs.append(job)

        if len(readset.annotation_files) > 0 and os.path.isfile(readset.annotation_files[0]):
            jobs.extend(self.verify_bam_id(readset, "DNA"))

        job = bvatools.depth_of_coverage(
            input,
            input_file_prefix + "metrics.targetCoverage.txt",
            full_coverage_bed,
            other_options=config.param('bvatools_depth_of_coverage', 'other_options', required=False),
            reference_genome=readset.reference_file
        )
        job.name = "bvatools_depth_of_coverage." + readset.name + ".doc" + "." + readset.run + "." + readset.lane
        job.samples = [readset.sample]
        job.report_files = [readset.bam + ".metrics.targetCoverage.txt"]
        jobs.append(job)

        return jobs

class RNARunProcessingAligner(RunProcessingAligner):
    def get_annotation_files(self):
        folder_name = os.path.basename(self.genome_folder)
        ini_file = os.path.join(self.genome_folder + os.sep + folder_name + ".ini")
        if os.path.isfile(ini_file):
            config.parse_files([ini_file])

            section = "DEFAULT"
            source = config.get(section, "source")
            version = config.get(section, "version")

            annotation_files = [
                os.path.join(self.genome_folder,
                             "annotations",
                             folder_name + '.' + source + version + ".transcript_id.gtf"),

                os.path.join(self.genome_folder,
                             "annotations",
                             "rrna_bwa_index",
                             folder_name + '.' + source + version + ".rrna.fa"),

                os.path.join(self.genome_folder,
                             "annotations",
                             folder_name + '.' + source + version + ".ref_flat.tsv")
            ]

            dbsnp_option_name = "dbsnp_version"
            af_option_name = "population_AF"
            if config.has_option(section, dbsnp_option_name) and config.has_option(section, af_option_name):
                dbsnp_version = config.get(section, dbsnp_option_name)
                af_name = config.get(section, af_option_name)
                annotation_files.append(
                    os.path.join(
                        self.genome_folder,
                        "annotations",
                        folder_name + ".dbSNP" + dbsnp_version + "_" + af_name + ".vcf.gz"
                    )
                )
            return annotation_files
        else:
            return None

    def get_metrics_jobs(self, readset):
        jobs = []
        jobs += self.verify_bam_id(readset, "RNA") + self._rnaseqc(readset) + self._picard_rna_metrics(readset) + \
                self._estimate_ribosomal_rna(readset, self.platform)
        return jobs

    @staticmethod
    def _rnaseqc(readset):
        """
        Computes a series of quality control metrics using
        [RNA-SeQC](https://www.broadinstitute.org/cancer/cga/rna-seqc).
        """
        jobs = []
        input_bam = readset.bam + ".dup.bam"
        input_bam_directory = os.path.dirname(input_bam)
        sample_file = input_bam + ".sample_file"
        sample_row = readset.sample.name + "\t" + input_bam + "\tRNAseq"
        output_directory = os.path.join(input_bam_directory, "rnaseqc_" + readset.sample.name + "." + readset.library)
        ribosomal_interval_file = os.path.join(output_directory, "empty.list")

        if len(readset.annotation_files) > 0 and os.path.isfile(readset.annotation_files[0]):
            gtf_transcript_id = readset.annotation_files[0]
            reference = readset.reference_file
            job = concat_jobs([
                bash.mkdir(output_directory),
                Job(command="touch " + ribosomal_interval_file),
                Job(
                    [input_bam],
                    [sample_file],
                    command="""\
echo "Sample\tBamFile\tNote\n{sample_row}" \\
 > {sample_file}""".format(
                        sample_row=sample_row,
                        sample_file=sample_file)),
                metrics.rnaseqc(
                    sample_file,
                    output_directory,
                    readset.fastq2 is None,
                    gtf_file=gtf_transcript_id,
                    ribosomal_interval_file=ribosomal_interval_file,
                    reference=reference
                ),
                bash.cp(
                    os.path.join(output_directory, "metrics.tsv"),
                    os.path.join(input_bam_directory, readset.sample.name + "." + readset.library + ".rnaseqc.sorted.dup.metrics.tsv")
                )],
                name="rnaseqc." + readset.name + ".rnaseqc" + "." + readset.run + "." + readset.lane,
                samples=[readset.sample],
                report_files=[os.path.join(input_bam_directory, readset.sample.name + "." + readset.library + ".rnaseqc.sorted.dup.metrics.tsv")]
            )
            jobs.append(job)

        return jobs

    @staticmethod
    def _picard_rna_metrics(readset):
        """
        Computes a series of quality control metrics using both CollectRnaSeqMetrics and CollectAlignmentSummaryMetrics
        functions metrics are collected using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        sample = readset.sample

        alignment_file = readset.bam + ".bam"
        output_directory = os.path.dirname(alignment_file)

        job = picard.collect_multiple_metrics(
            alignment_file, readset.bam + ".metrics",
            reference_sequence=readset.reference_file
        )
        job.name = "picard_collect_multiple_metrics." + readset.name + ".met" + "." + readset.run + "." + readset.lane
        job.samples = [readset.sample]
        job.report_files = [
            readset.bam + ".metrics.alignment_summary_metrics",
            readset.bam + ".metrics.insert_size_metrics"
        ]
        jobs.append(job)

        if len(readset.annotation_files) > 2 and os.path.isfile(readset.annotation_files[2]):
            job = picard.collect_rna_metrics(
                alignment_file,
                os.path.join(output_directory, sample.name),
                readset.annotation_files[2],
                readset.reference_file
            )

            job.name = "picard_rna_metrics." + readset.name + ".rmet" + "." + readset.run + "." + readset.lane
            job.samples = [readset.sample]
            jobs.append(job)

        return jobs

    @staticmethod
    def _estimate_ribosomal_rna(readset, platform):
        """
        Use bwa mem to align reads on the rRNA reference fasta and count the number of read mapped
        The filtered reads are aligned to a reference fasta file of ribosomal sequence. The alignment is done per
        sequencing readset.
        The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
        BWA output BAM files are then sorted by coordinate using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        if len(readset.annotation_files) > 1 and os.path.isfile(readset.annotation_files[0]) and os.path.isfile(
                readset.annotation_files[1]):
            readset_bam = readset.bam + ".bam"
            readset_metrics_bam = readset.bam + ".rRNA.bam"

            job = concat_jobs([
                pipe_jobs([
                    bvatools.bam2fq(readset_bam),
                    bwa.mem(
                        "/dev/stdin",
                        None,
                        read_group=RunProcessingAligner.get_rg_tag(readset, platform, 'bwa_mem_rRNA'),
                        ref=readset.annotation_files[1],
                        ini_section='bwa_mem_rRNA'
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        readset_metrics_bam,
                        "coordinate",
                        ini_section='picard_sort_sam_rrna'
                    )
                ]),
                tools.py_rrnaBAMcount(
                    bam=readset_metrics_bam,
                    gtf=readset.annotation_files[0],
                    output=os.path.join(readset.bam + ".metrics.rRNA.tsv"),
                    typ="transcript"
                )],
                name="bwa_mem_rRNA." + readset.name + ".rRNA" + "." + readset.run + "." + readset.lane,
                samples=[readset.sample],
                report_files=[os.path.join(readset.bam + ".metrics.rRNA.tsv")],
                removable_files=[readset_metrics_bam]
            )

            jobs.append(job)

        return jobs

class StarRunProcessingAligner(RNARunProcessingAligner):
    def __init__(self, output_dir, genome_folder, nb_cycles, platform):
        super(StarRunProcessingAligner, self).__init__(output_dir, genome_folder, platform)
        self._nb_cycles = nb_cycles

    @property
    def nb_cycles(self):
        return self._nb_cycles

    def get_reference_index(self):
        folder_name = os.path.basename(self.genome_folder)
        ini_file = os.path.join(self.genome_folder + os.sep + folder_name + ".ini")
        if os.path.isfile(ini_file):
            config.parse_files([ini_file])

            source = config.get("DEFAULT", "source")
            version = config.get("DEFAULT", "version")
            star_version=config.get("DEFAULT", "module_star").split('/')[-1]

            star_index_folder = os.path.join(self.genome_folder, "genome", "star" + star_version + "_index")
            if not os.path.exists(os.path.expandvars(star_index_folder)):
                star_index_folder = os.path.join(self.genome_folder, "genome", "star_index")
                if not os.path.exists(os.path.expandvars(star_index_folder)):
                    return None

            indexFile = os.path.join(star_index_folder, source + version + ".sjdbOverhang" + str(self.nb_cycles - 1))
            indexFile2 = os.path.join(star_index_folder, source + version + ".sjdbOverhang" + str(self.nb_cycles - 2))
            return indexFile if os.path.exists(indexFile) else indexFile2
        else:
            return None

    def get_reference_file(self):
        folder_name = os.path.basename(self.genome_folder)
        return os.path.join(
            self.genome_folder,
            "genome",
            folder_name + ".fa"
        )

    def get_alignment_job(self, readset):
        output = readset.bam + ".bam"

        rg_center = config.param('star_align', 'sequencing_center', required=False)

        # We can't set the exact bam filename for STAR, so we output the result in a specific directory, and move the
        # bam to the expected place with the right name.
        star_bam_name = "Aligned.sortedByCoord.out.bam"
        star_output_directory = os.path.join(os.path.dirname(output), readset.library)

        star_job = star.align(
            reads1=readset.fastq1,
            reads2=readset.fastq2,
            output_directory=star_output_directory,
            sort_bam=True,
            genome_index_folder=readset.aligner_reference_index,
            rg_id=readset.library + "_" + readset.run + "_" + readset.lane,
            rg_sample=readset.sample.name,
            rg_library=readset.library if readset.library else "",
            rg_platform_unit=readset.run + "_" + readset.lane if readset.run and readset.lane else "",
            rg_platform=self.platform,
            rg_center=rg_center if rg_center else ""
        )
        # we clean the output of the star job since we move the file, and the moved file is the output of the move job
        star_job._output_files = []

        move_job = bash.mv(
            os.path.join(star_output_directory, star_bam_name),
            output
        )
        # we clean the input of the move job to be conistent with te star job
        move_job._input_files = []

        job = concat_jobs(
            [
                star_job,
                move_job,
                picard.build_bam_index(
                    output,
                    output[::-1].replace(".bam"[::-1], ".bai"[::-1], 1)[::-1]
                )
            ],
            name="star_align." + readset.name + "." + readset.run + "." + readset.lane,
            samples=[readset.sample]
        )
        return job


class CellrangerRunProcessingAligner(RNARunProcessingAligner):
    def get_reference_index(self):
        assembly_name = os.path.basename(self.genome_folder)
        ini_file = os.path.join(self.genome_folder, assembly_name + ".ini")
        if os.path.isfile(ini_file):
            config.parse_files([ini_file])

            if config.has_option("DEFAULT", "10x_transcriptome"):
                transcriptome = config.get("DEFAULT", "10x_transcriptome")
                trans_path = os.path.join(self.genome_folder, "genome", "10xGenomics", transcriptome)

                if trans_path is None or not os.path.isdir(trans_path):
                    """ Main transcriptome not found, we will try to resolve it if synomym exists """
                    if config.has_option("DEFAULT", "assembly_synonyms"):
                        synonym = config.get("DEFAULT", "assembly_synonyms")
                        trans_path = self.find_10x_synonym_reference(self.genome_folder, synonym, "10x_transcriptome")

            return trans_path
        else:
            return ""

    def get_reference_file(self):
        return os.path.join(self.get_reference_index(), "fasta", "genome.fa")

    def get_alignment_job(self, readset):
        cell_count_job = cellranger.count(
            read1=readset.fastq1,
            read2=readset.fastq2,
            sample_id=readset_name + "_" + readset.sample_number,
            lane_id=readset.lane,
            project=readset.project,
            fastqs=os.path.dirname(readset.fastq1),
            transcriptome=self.get_reference_index()
        )

        outdir = os.path.join(
            self.output_dir,
            config.param('cellranger_count', 'working_dir', required=True) + "_" + readset.lane_str,
            readset.name + "_" + readset.sample_number,
            "outs"
        )
        default_bam_name = "possorted_genome_bam"
        input_bam = os.path.join(outdir, default_bam_name) + ".bam"
        input_bai = os.path.join(outdir, default_bam_name) + ".bam.bai"
        input_summary = os.path.join(outdir, "web_summary.html")
        output_bam = readset.bam + ".bam"
        output_bai = readset.bam + ".bai"
        output_summary = readset.bam + ".10x_summary.html"
        output_zip = readset.bam + "." + readset.run + "." + readset.lane_str + ".10x_outputs.zip"
        job = concat_jobs(
            [
                bash.mkdir(os.path.dirname(readset.bam + ".bam")),
                cell_count_job,
                bash.cp(
                    input_bam,
                    output_bam
                ),
                bash.cp(
                    input_bai,
                    output_bai
                ),
                bash.cp(
                    input_summary,
                    output_summary
                ),
                Job(
                    output_files=[output_zip],
                    command="cd " + outdir + " && zip -r " + output_zip + " *"
                )
            ],
            name="cellranger_count." + readset.name + "." + readset.run + "." + readset.lane,
            samples=[readset.sample]
        )
        return job

class AtacRunProcessingAligner(RNARunProcessingAligner):
    def get_reference_index(self):
        assembly_name = os.path.basename(self.genome_folder)
        ini_file = os.path.join(self.genome_folder + os.sep + assembly_name + ".ini")
        if os.path.isfile(ini_file):
            config.parse_files([ini_file])
            trans_path = ""

            if config.has_option("DEFAULT", "10x_atac_transcriptome"):
                transcriptome = config.get("DEFAULT", "10x_atac_transcriptome")
                trans_path = os.path.join(self.genome_folder, "genome", "10xGenomics", transcriptome)

            if trans_path is None or not os.path.isdir(trans_path):
                """ Main transcriptome not found, we will try to resolve it if synomym exists """
                if config.has_option("DEFAULT", "assembly_synonyms"):
                    synonym = config.get("DEFAULT", "assembly_synonyms")
                    trans_path = self.find_10x_synonym_reference(self.genome_folder, synonym, "10x_atac_transcriptome")
            return trans_path
        else:
            return ""

    def get_reference_file(self):
        return os.path.join(self.get_reference_index(), "fasta", "genome.fa")

    def get_alignment_job(self, readset):
        cell_atac_job = cellranger.atac(
            reads1=readset.fastq1,
            reads2=readset.fastq2,
            sample_id=readset.name + "_" + readset.sample_number,
            lane_id=readset.lane,
            project=readset.project,
            fastqs=os.path.dirname(readset.fastq1),
            reference=self.get_reference_index()
        )
        outdir = os.path.join(
            self.output_dir,
            config.param('cellranger_atac', 'working_dir', required=True) + "_" + readset.lane,
            readset.name + "_" + readset.sample_number,
            "outs"
        )
        default_bam_name = "possorted_bam"
        input_bam = os.path.join(outdir, default_bam_name) + ".bam"
        input_bai = os.path.join(outdir, default_bam_name) + ".bam.bai"
        input_summary = os.path.join(outdir, "web_summary.html")
        output_bam = readset.bam + ".bam"
        output_bai = readset.bam + ".bai"
        output_summary = readset.bam + ".10x_summary.html"
        output_zip = readset.bam + "." + readset.run + "." + readset.lane + ".10x_outputs.zip"
        job = concat_jobs(
            [
                bash.mkdir(os.path.dirname(readset.bam + ".bam")),
                cell_atac_job,
                bash.cp(
                    input_bam,
                    output_bam
                ),
                bash.cp(
                    input_bai,
                    output_bai
                ),
                bash.cp(
                    input_summary,
                    output_summary
                ),
                Job(
                    output_files=[output_zip],
                    command="cd " + outdir + " && zip -r " + output_zip + " *"
                )
            ],
            name="cellranger_atac." + readset.name + "." + readset.run + "." + readset.lane,
            samples=[readset.sample]
        )
        return job

class CellCounterRunProcessingAligner(NullRunProcessingAligner):
    """ Fake aligner used when a 10x single cell run has only a single read. We only count UMIs by barcode. """

    def get_fastq_metrics_jobs(self, readset):
        input_file = readset.fastq1
        output = os.path.join(self.output_dir,
                              "Unaligned." + str(readset.lane_str) + "Count",
                              readset.name + '_S' + readset.sample_number + '_L00' + readset.lane_str
                              + '_R1_001.count.csv')
        job = Job([input_file], [output], [["cell_counter", "module_java"]],
                  name="cell_counter." + readset.name + "." + readset.run + "." + readset.lane_str)
        job.command = """\java {java_other_options} -cp {jar} ca.mcgill.genome.mps.core.util.ChromiumCellCounter {input} {output}""".format(
            java_other_options=config.param('cell_counter', 'java_other_options'),
            jar=config.param('cell_counter', 'jar'),
            input=input_file,
            output=output
        )

        return [job]

