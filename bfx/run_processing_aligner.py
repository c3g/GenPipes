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

import configparser
import os
import re

from core.job import concat_jobs, pipe_jobs, Job
from core.config import global_config_parser
from bfx import bvatools
from bfx import snpeff
from bfx import verify_bam_id
from bfx import bwa
from bfx import metrics
from bfx import picard
from bfx import star
from bfx import tools


class RunProcessingAligner(object):
    def __init__(self, output_dir, genome_folder):
        self._output_dir = output_dir
        self._genome_folder = genome_folder

    @property
    def output_dir(self):
        return self._output_dir

    @property
    def genome_folder(self):
        return self._genome_folder

    def get_reference_index(self):
        raise NotImplementedError("Please Implement this method")

    def get_alignment_job(self, readset):
        raise NotImplementedError("Please Implement this method")

    def get_metrics_jobs(self, readset):
        raise NotImplementedError("Please Implement this method")

    def get_annotation_files(self):
        raise NotImplementedError("Please Implement this method")

    @staticmethod
    def get_rg_tag(readset, ini_section):
        return "'@RG" + \
               "\tID:" + readset.library + "_" + readset.run + "_" + readset.lane + \
               "\tSM:" + readset.sample.name + \
               "\tLB:" + readset.library + \
               "\tPU:run" + readset.run + "_" + readset.lane + \
               ("\tCN:" + global_config_parser.param(ini_section, 'sequencing_center')
                if global_config_parser.param(ini_section, 'sequencing_center', required=False) else "") + \
               "\tPL:Illumina" + \
               "'"


class BwaRunProcessingAligner(RunProcessingAligner):
    downloaded_bed_files = []
    created_interval_lists = []

    def get_reference_index(self):
        folder_name = os.path.basename(self.genome_folder)
        return os.path.join(self.genome_folder,
                            "genome",
                            "bwa_index",
                            folder_name + ".fa")

    def get_annotation_files(self):
        folder_name = os.path.basename(self.genome_folder)
        ini_file = os.path.join(self.genome_folder + os.sep + folder_name + ".ini")
        if os.path.isfile(ini_file):
            genome_config = configparser.SafeConfigParser()
            genome_config.read(ini_file)

            section = "DEFAULT"
            dbsnp_option_name = "dbsnp_version"
            af_option_name = "population_AF"

            if genome_config.has_option(section, dbsnp_option_name) and\
                    genome_config.has_option(section, af_option_name):
                dbsnp_version = genome_config.get(section, dbsnp_option_name)
                af_name = genome_config.get(section, af_option_name)
                return [
                    os.path.join(self.genome_folder,
                                 "annotations",
                                 folder_name + ".dbSNP" + dbsnp_version + "_" + af_name + ".vcf"),
                ]

        return []

    def get_alignment_job(self, readset):
        output = readset.bam + ".bam"
        job = concat_jobs([
                              Job(command="mkdir -p " + os.path.dirname(output)),
                              pipe_jobs([
                                  bwa.mem(
                                      readset.fastq1,
                                      readset.fastq2,
                                      read_group=RunProcessingAligner.get_rg_tag(readset, 'bwa_mem'),
                                      ref=readset.aligner_reference_index
                                  ),
                                  picard.sort_sam(
                                      "/dev/stdin",
                                      output,
                                      "coordinate"
                                  )
                              ])
                          ], name="bwa_mem_picard_sort_sam." + readset.name + "." + readset.run + "." + readset.lane)

        return job

    def get_metrics_jobs(self, readset):
        jobs = []

        input_file_prefix = readset.bam + '.'
        input = input_file_prefix + "bam"

        job = picard.collect_multiple_metrics(input, input_file_prefix + "metrics",
                                              reference_sequence=readset.reference_file)
        job.name = "picard_collect_multiple_metrics." + readset.name + ".met" + "." + readset.run + "." + readset.lane
        jobs.append(job)

        if readset.beds:
            coverage_bed = readset.beds[0]
            full_coverage_bed = (self.output_dir + os.sep + coverage_bed)
        else:
            coverage_bed = None
            full_coverage_bed = None

        if coverage_bed:
            if (not os.path.exists(full_coverage_bed)) and \
                    (coverage_bed not in BwaRunProcessingAligner.downloaded_bed_files):
                # Download the bed file
                command = global_config_parser.param('DEFAULT', 'fetch_bed_file_command').format(
                    output_directory=self.output_dir,
                    filename=coverage_bed
                )
                job = Job([], [full_coverage_bed], command=command, name="bed_download." + coverage_bed)
                BwaRunProcessingAligner.downloaded_bed_files.append(coverage_bed)
                jobs.append(job)

            interval_list = re.sub("\.[^.]+$", ".interval_list", coverage_bed)

            if interval_list not in BwaRunProcessingAligner.created_interval_lists:
                # Create one job to generate the interval list from the bed file
                ref_dict = os.path.splitext(readset.reference_file)[0] + '.dict'
                job = tools.bed2interval_list(ref_dict, full_coverage_bed, interval_list)
                job.name = "interval_list." + coverage_bed
                BwaRunProcessingAligner.created_interval_lists.append(interval_list)
                jobs.append(job)

            job = picard.calculate_hs_metrics(input_file_prefix + "bam", input_file_prefix + "metrics.onTarget.txt",
                                              interval_list, reference_sequence=readset.reference_file)
            job.name = "picard_calculate_hs_metrics." + readset.name + ".hs" + "." + readset.run + "." + readset.lane
            jobs.append(job)

        jobs.extend(self.verify_bam_id(readset))

        job = bvatools.depth_of_coverage(
            input,
            input_file_prefix + "metrics.targetCoverage.txt",
            full_coverage_bed,
            other_options=global_config_parser.param('bvatools_depth_of_coverage', 'other_options', required=False),
            reference_genome=readset.reference_file
        )
        job.name = "bvatools_depth_of_coverage." + readset.name + ".doc" + "." + readset.run + "." + readset.lane
        jobs.append(job)

        return jobs

    def verify_bam_id(self, readset):
        """
            verifyBamID is a software that verifies whether the reads in particular file match previously known
            genotypes for an individual (or group of individuals), and checks whether the reads are contaminated
            as a mixture of two samples. verifyBamID can detect sample contamination and swaps when external
            genotypes are available. When external genotypes are not available, verifyBamID still robustly
            detects sample swaps.
        """
        jobs = []
        if len(readset.annotation_files) > 0 and os.path.isfile(readset.annotation_files[0]):

            known_variants_annotated = readset.annotation_files[0]
            known_variants_annotated_filtered = known_variants_annotated

            input_bam = readset.bam + ".bam"
            output_prefix = readset.bam + ".metrics.verifyBamId"

            jobs.append(concat_jobs([
                verify_bam_id.verify(
                    input_bam,
                    known_variants_annotated_filtered,
                    output_prefix,
                ),
                # the first column starts with a # (a comment for nanuq) so we remove the column and output the result
                # in a file with the name supported by nanuq
                Job([output_prefix + ".selfSM"],
                    [output_prefix + ".tsv"],
                    command="cut -f2- " + output_prefix + ".selfSM > " + output_prefix + ".tsv"
                )
                ],
                name="verify_bam_id." + readset.name + "." + readset.run + "." + readset.lane
            ))

        return jobs


class StarRunProcessingAligner(RunProcessingAligner):
    def __init__(self, output_dir, genome_folder, nb_cycles):
        super(StarRunProcessingAligner, self).__init__(output_dir, genome_folder)
        self._nb_cycles = nb_cycles

    @property
    def nb_cycles(self):
        return self._nb_cycles

    def get_reference_index(self):
        folder_name = os.path.basename(self.genome_folder)
        ini_file = os.path.join(self.genome_folder + os.sep + folder_name + ".ini")
        if os.path.isfile(ini_file):
            genome_config = configparser.SafeConfigParser()
            genome_config.read(ini_file)

            source = genome_config.get("DEFAULT", "source")
            version = genome_config.get("DEFAULT", "version")

            return os.path.join(self.genome_folder,
                                "genome",
                                "star_index",
                                source + version + ".sjdbOverhang" + str(self.nb_cycles - 1))
        else:
            return None

    def get_annotation_files(self):
        folder_name = os.path.basename(self.genome_folder)
        ini_file = os.path.join(self.genome_folder + os.sep + folder_name + ".ini")
        if os.path.isfile(ini_file):
            genome_config = configparser.SafeConfigParser()
            genome_config.read(ini_file)

            source = genome_config.get("DEFAULT", "source")
            version = genome_config.get("DEFAULT", "version")

            return [
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
        else:
            return None

    def get_alignment_job(self, readset):
        output = readset.bam + ".bam"

        rg_center = global_config_parser.param('star_align', 'sequencing_center', required=False)

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
            rg_platform="Illumina",
            rg_center=rg_center if rg_center else ""
        )
        # we clean the output of the star job since we move the file, and the moved file is the output of the move job
        star_job._output_files = []

        job = concat_jobs([
            star_job,
            Job(output_files=[output], command="mv " + os.path.join(star_output_directory, star_bam_name) + " "
                                               + output),
            picard.build_bam_index(output, output[::-1].replace(".bam"[::-1], ".bai"[::-1], 1)[::-1])
        ])
        job.name = "star_align." + readset.name + "." + readset.run + "." + readset.lane
        return job

    def get_metrics_jobs(self, readset):
        jobs = []
        jobs += StarRunProcessingAligner._rnaseqc(readset) + StarRunProcessingAligner._picard_rna_metrics(readset) + \
                StarRunProcessingAligner._estimate_ribosomal_rna(readset)
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
                                  Job(command="mkdir -p " + output_directory),
                                  Job(command="touch " + ribosomal_interval_file),
                                  Job([input_bam], [sample_file], command="""\
        echo "Sample\tBamFile\tNote
        {sample_row}" \\
        > {sample_file}""".format(sample_row=sample_row, sample_file=sample_file)),
                                  metrics.rnaseqc(sample_file, output_directory, readset.fastq2 is not None,
                                                  gtf_file=gtf_transcript_id,
                                                  ribosomal_interval_file=ribosomal_interval_file, reference=reference),
                                  Job(command="cp " + os.path.join(output_directory,
                                                                   "metrics.tsv") + " " + os.path.join(
                                      input_bam_directory,
                                      readset.sample.name + "." +
                                      readset.library +
                                      ".rnaseqc.sorted.dup.metrics.tsv"))
                              ], name="rnaseqc." + readset.name + ".rnaseqc" + "." + readset.run + "." + readset.lane)
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

        job = picard.collect_multiple_metrics(alignment_file, readset.bam + ".metrics",
                                              reference_sequence=readset.reference_file)
        job.name = "picard_collect_multiple_metrics." + readset.name + ".met" + "." + readset.run + "." + readset.lane
        jobs.append(job)

        if len(readset.annotation_files) > 2 and os.path.isfile(readset.annotation_files[2]):
            job = picard.collect_rna_metrics(alignment_file,
                                             os.path.join(output_directory, sample.name),
                                             readset.annotation_files[2],
                                             readset.reference_file)

            job.name = "picard_rna_metrics." + readset.name + ".rmet" + "." + readset.run + "." + readset.lane
            jobs.append(job)

        return jobs

    @staticmethod
    def _estimate_ribosomal_rna(readset):
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
                                      bvatools.bam2fq(
                                          readset_bam
                                      ),
                                      bwa.mem(
                                          "/dev/stdin",
                                          None,
                                          read_group=RunProcessingAligner.get_rg_tag(readset, 'bwa_mem_rRNA'),
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
                                      typ="transcript")],
                              name="bwa_mem_rRNA." + readset.name + ".rRNA" + "." + readset.run + "." + readset.lane)

            job.removable_files = [readset_metrics_bam]
            jobs.append(job)

        return jobs


