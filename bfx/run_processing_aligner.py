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

from core.job import *
from core.config import *
from bfx import bvatools
from bfx import bwa
from bfx import metrics
from bfx import picard
from bfx import star
from bfx import tools


class RunProcessingAligner(object):
    def __init__(self, output_dir):
        self._output_dir = output_dir

    @property
    def output_dir(self):
        return self._output_dir

    def get_reference_index(self, genome_folder):
        raise NotImplementedError("Please Implement this method")

    def get_alignment_jobs(self, readset):
        raise NotImplementedError("Please Implement this method")

    def get_metrics_jobs(self, readset):
        raise NotImplementedError("Please Implement this method")

    def get_annotation_files(self, genome_folder):
        raise NotImplementedError("Please Implement this method")

    @staticmethod
    def get_rg_tag(readset, ini_section):
        return "'@RG" + \
               "\tID:" + readset.library + "_" + readset.run + "_" + readset.lane + \
               "\tSM:" + readset.sample.name + \
               "\tLB:" + readset.library + \
               "\tPU:run" + readset.run + "_" + readset.lane + \
               ("\tCN:" + config.param(ini_section, 'sequencing_center')
                if config.param(ini_section, 'sequencing_center', required=False) else "") + \
               "\tPL:Illumina" + \
               "'"


class BwaRunProcessingAligner(RunProcessingAligner):
    @property
    def created_interval_lists(self):
        if not hasattr(self, "_created_interval_lists"):
            self._created_interval_lists = []
        return self._created_interval_lists

    @property
    def downloaded_bed_files(self):
        if not hasattr(self, "_downloaded_bed_files"):
            self._downloaded_bed_files = []
        return self._downloaded_bed_files

    def get_reference_index(self, genome_folder):
        folder_name = os.path.basename(genome_folder)
        return os.path.join(genome_folder,
                            "genome",
                            "bwa_index",
                            folder_name + ".fa")

    def get_annotation_files(self, genome_folder):
        return []

    def get_alignment_jobs(self, readset):
        jobs = []
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
                          ], name="bwa_mem_picard_sort_sam." + readset.name + "_" + readset.run + "_" + readset.lane)

        jobs.append(job)
        return jobs

    def get_metrics_jobs(self, readset):
        jobs = []

        input_file_prefix = readset.bam + '.'
        input = input_file_prefix + "bam"

        job = picard.collect_multiple_metrics(input, input_file_prefix + "metrics",
                                              reference_sequence=readset.reference_file)
        job.name = "picard_collect_multiple_metrics." + readset.name + ".met" + "." + readset.run + "." + readset.lane
        jobs.append(job)

        coverage_bed = bvatools.resolve_readset_coverage_bed(readset)
        full_coverage_bed = (self.output_dir + os.sep + coverage_bed) if coverage_bed else None

        if coverage_bed:
            if (not os.path.exists(full_coverage_bed)) and (coverage_bed not in self.downloaded_bed_files):
                # Download the bed file
                command = config.param('DEFAULT', 'fetch_bed_file_command').format(
                    output_directory=self.output_dir,
                    filename=coverage_bed
                )
                job = Job([], [full_coverage_bed], command=command, name="bed_download." + coverage_bed)
                self.downloaded_bed_files.append(coverage_bed)
                jobs.append(job)

            interval_list = re.sub("\.[^.]+$", ".interval_list", coverage_bed)

            if interval_list not in self.created_interval_lists:
                # Create one job to generate the interval list from the bed file
                ref_dict = os.path.splitext(readset.reference_file)[0] + '.dict'
                job = tools.bed2interval_list(ref_dict, full_coverage_bed, interval_list)
                job.name = "interval_list." + coverage_bed
                self.created_interval_lists.append(interval_list)
                jobs.append(job)

            job = picard.calculate_hs_metrics(input_file_prefix + "bam", input_file_prefix + "metrics.onTarget.txt",
                                              interval_list, reference_sequence=readset.reference_file)
            job.name = "picard_calculate_hs_metrics." + readset.name + ".hs" + "." + readset.run + "." + readset.lane
            jobs.append(job)

        job = bvatools.depth_of_coverage(
            input,
            input_file_prefix + "metrics.targetCoverage.txt",
            full_coverage_bed,
            other_options=config.param('bvatools_depth_of_coverage', 'other_options', required=False),
            reference_genome=readset.reference_file
        )
        job.name = "bvatools_depth_of_coverage." + readset.name + ".doc" + "." + readset.run + "." + readset.lane
        jobs.append(job)

        return jobs


class StarRunProcessingAligner(RunProcessingAligner):
    def __init__(self, output_dir, nb_cycles):
        super(StarRunProcessingAligner, self).__init__(output_dir)
        self._nb_cycles = nb_cycles

    @property
    def nb_cycles(self):
        return self._nb_cycles

    def get_reference_index(self, genome_folder):
        folder_name = os.path.basename(genome_folder)
        ini_file = os.path.join(genome_folder + os.sep + folder_name + ".ini")
        genome_config = ConfigParser.SafeConfigParser()
        genome_config.read(ini_file)

        source = genome_config.get("DEFAULT", "source")
        version = genome_config.get("DEFAULT", "version")

        return os.path.join(genome_folder,
                            "genome",
                            "star_index",
                            source + version + ".sjdbOverhang" + str(self.nb_cycles - 1))

    def get_annotation_files(self, genome_folder):
        folder_name = os.path.basename(genome_folder)
        ini_file = os.path.join(genome_folder + os.sep + folder_name + ".ini")
        genome_config = ConfigParser.SafeConfigParser()
        genome_config.read(ini_file)

        source = genome_config.get("DEFAULT", "source")
        version = genome_config.get("DEFAULT", "version")

        return [
            os.path.join(genome_folder,
                         "annotations",
                         folder_name + '.' + source + version + ".transcript_id.gtf"),

            os.path.join(genome_folder,
                         "annotations",
                         "rrna_bwa_index",
                         folder_name + '.' + source + version + ".rrna.fa"),

            os.path.join(genome_folder,
                         "annotations",
                         folder_name + '.' + source + version + ".ref_flat.tsv")
        ]

    def get_alignment_jobs(self, readset):
        jobs = []
        output = readset.bam + ".bam"

        rg_center = config.param('star_align', 'sequencing_center', required=False)
        star_bam_name = "Aligned.sortedByCoord.out.bam"

        job = concat_jobs([
            star.align(
                reads1=readset.fastq1,
                reads2=readset.fastq2,
                output_directory=os.path.dirname(output),
                sort_bam=True,
                genome_index_folder=readset.aligner_reference_index,
                rg_id=readset.library + "_" + readset.run + "_" + readset.lane,
                rg_sample=readset.sample.name,
                rg_library=readset.library if readset.library else "",
                rg_platform_unit=readset.run + "_" + readset.lane if readset.run and readset.lane else "",
                rg_platform="Illumina",
                rg_center=rg_center if rg_center else ""
            ),
            Job(output_files=[output], command="mv " + os.path.dirname(output) + os.sep + star_bam_name + " " + output),
            Job(command="ln -s " + output + " " + os.path.dirname(output) + os.sep + star_bam_name),
            picard.build_bam_index(output, output[::-1].replace(".bam"[::-1], ".bai"[::-1], 1)[::-1])
        ])
        job.name = "star_align." + readset.name + "." + readset.run + "." + readset.lane
        jobs.append(job)
        return jobs

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
                              ], name="rnaseqc" + readset.name + ".rnaseqc" + "." + readset.run + "." + readset.lane)
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


