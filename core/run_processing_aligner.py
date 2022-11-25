################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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
import configparser

# MUGQIC Modules
from .job import Job, concat_jobs, pipe_jobs
from .config import config, _raise, SanitycheckError
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

    def get_fastq_metrics_jobs(self, readset):
        return[]

    def get_annotation_files(self):
        raise NotImplementedError("Please Implement this method")

    def is_run_mark_duplicate(self):
        return True

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

    def verify_bam_id(self, readset):
        """
            verifyBamID is a software that verifies whether the reads in particular file match previously known
            genotypes for an individual (or group of individuals), and checks whether the reads are contaminated
            as a mixture of two samples. verifyBamID can detect sample contamination and swaps when external
            genotypes are available. When external genotypes are not available, verifyBamID still robustly
            detects sample swaps.
        """
        jobs = []

        known_variants_annotated_filtered = ""
        if verify_bam_id.getVersion() == 1:
            if 'vcf' in readset.annotation_files and readset.annotation_files['vcf'] and os.path.isfile(readset.annotation_files['vcf']):
                known_variants_annotated_filtered = readset.annotation_files['vcf']
            else:
                log.info("VerifyBamID is requested but no suitable VCF was found for " + readset.species + "... skipping VerifyBamID...")
        else:
            if 'dat' in readset.annotation_files and os.path.isfile(readset.annotation_files['dat']+".UD") and os.path.isfile(readset.annotation_files['dat']+".bed") and os.path.isfile(readset.annotation_files['dat']+".mu"):
                known_variants_annotated_filtered = readset.annotation_files['dat']
            else:
                log.info("VerifyBamID2 is requested but no SVD dataset was found for " + readset.species + "... skipping VerifyBamID...")

        if known_variants_annotated_filtered:
            input_bam = readset.bam + ".bam"
            output_prefix = readset.bam + ".metrics.verifyBamId"
            jobs.append(
                concat_jobs(
                    [
                        verify_bam_id.verify(
                            input_bam,
                            output_prefix,
                            var=known_variants_annotated_filtered,
                            ref=readset.reference_file
                        )
                    ],
                    name="verify_bam_id." + readset.name + "." + readset.run + "." + readset.lane,
                    report_files=[readset.bam + ".metrics.verifyBamId.selfSM"],
                    samples=[readset.sample]
                )
            )

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
        return {}

class BwaRunProcessingAligner(RunProcessingAligner):
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

    def get_dictionary_file(self):
        folder_name = os.path.basename(self.genome_folder)
        return os.path.join(
            self.genome_folder,
            "genome",
            folder_name + ".dict"
        )

    def get_annotation_files(self):
        folder_name = os.path.basename(self.genome_folder)
        ini_file = os.path.join(self.genome_folder + os.sep + folder_name + ".ini")
        if os.path.isfile(ini_file):
            # config.parse_files([ini_file])
            genome_config = configparser.SafeConfigParser()
            genome_config.read(ini_file)

            section = "DEFAULT"
            dbsnp_option_name = "dbsnp_version"
            af_option_name = "population_AF"

            annotation_files = {}
            if genome_config.has_option(section, dbsnp_option_name) and genome_config.has_option(section, af_option_name):
                dbsnp_version = genome_config.get(section, dbsnp_option_name)
                af_name = genome_config.get(section, af_option_name)
                annotation_folder = os.path.join(self.genome_folder, "annotations")
                annotation_files["vcf"] = os.path.join(annotation_folder, folder_name + ".dbSNP" + dbsnp_version + "_" + af_name + ".vcf.gz")

                if verify_bam_id.getVersion() == 2:
                    annotation_files["dat"] = os.path.join(annotation_folder, "svd_datasets", folder_name + ".1000g.phase3.100k.vcf.gz.dat")
                    # log.debug(annotation_files["dat"])

        return annotation_files

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
            if (not os.path.exists(full_coverage_bed)):
                _raise(SanitycheckError("Could not find the bed coverage file " + full_coverage_bed + "..."))

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

        jobs.extend(self.verify_bam_id(readset))

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

            annotation_folder = os.path.join(self.genome_folder, "annotations")

            annotation_files = {
                "gtf": os.path.join(annotation_folder, folder_name + '.' + source + version + ".transcript_id.gtf"),
                "rrna_fa": os.path.join(annotation_folder, "rrna_bwa_index", folder_name + '.' + source + version + ".rrna.fa"),
                "ref_flat": os.path.join(annotation_folder, folder_name + '.' + source + version + ".ref_flat.tsv")
            }

            dbsnp_option_name = "dbsnp_version"
            af_option_name = "population_AF"
            if config.has_option(section, dbsnp_option_name) and config.has_option(section, af_option_name):
                dbsnp_version = config.get(section, dbsnp_option_name)
                af_name = "_" + config.get(section, af_option_name)
                annotation_files['vcf'] = os.path.join(annotation_folder, folder_name + ".dbSNP" + dbsnp_version + af_name + ".vcf.gz")

                if verify_bam_id.getVersion() == 2:
                    annotation_files["dat"] = os.path.join(annotation_folder, "svd_datasets", folder_name + ".1000g.phase3.100k.vcf.gz.dat")

            return annotation_files

    def get_metrics_jobs(self, readset):
        jobs = []
        jobs += self.verify_bam_id(readset) + self._rnaseqc(readset) + self._picard_rna_metrics(readset) + \
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

        if readset.annotation_files['gtf'] and os.path.isfile(readset.annotation_files['gtf']):
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
                    gtf_file=readset.annotation_files['gtf'],
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
            alignment_file,
            readset.bam + ".metrics",
            reference_sequence=readset.reference_file
        )
        job.name = "picard_collect_multiple_metrics." + readset.name + ".met" + "." + readset.run + "." + readset.lane
        job.samples = [readset.sample]
        job.report_files = [
            readset.bam + ".metrics.alignment_summary_metrics",
            readset.bam + ".metrics.insert_size_metrics"
        ]
        jobs.append(job)

        if readset.annotation_files['ref_flat'] and os.path.isfile(readset.annotation_files['ref_flat']):
            job = picard.collect_rna_metrics(
                alignment_file,
                os.path.join(output_directory, sample.name),
                readset.annotation_files['ref_flat'],
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
        if readset.annotation_files['gtf'] and readset.annotation_files['rrna_fa'] and os.path.isfile(readset.annotation_files['gtf']) and os.path.isfile(readset.annotation_files['rrna_fa']):
            readset_bam = readset.bam + ".bam"
            readset_metrics_bam = readset.bam + ".rRNA.bam"

            job = concat_jobs([
                pipe_jobs([
                    bvatools.bam2fq(readset_bam),
                    bwa.mem(
                        "/dev/stdin",
                        None,
                        read_group=RunProcessingAligner.get_rg_tag(readset, platform, 'bwa_mem_rRNA'),
                        ref=readset.annotation_files['rrna_fa'],
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
                    gtf=readset.annotation_files['gtf'],
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

    def get_dictionary_file(self):
        folder_name = os.path.basename(self.genome_folder)
        return os.path.join(
            self.genome_folder,
            "genome",
            folder_name + ".dict"
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

    def get_dictionary_file(self):
        folder_name = os.path.basename(self.genome_folder)
        return os.path.join(
            self.genome_folder,
            "genome",
            folder_name + ".dict"
        )

    def get_alignment_job(self, readset):
        outdir = os.path.join(
            self.output_dir,
            "10x_cellcount." + readset.lane,
            readset.name + "_" + readset.sample_number,
            "outs"
        )
        default_bam_name = "possorted_genome_bam"
        tmp_bam = os.path.join(outdir, readset.name + "." + default_bam_name + ".bam")
        tmp_bai = os.path.join(outdir, readset.name + "." + default_bam_name + ".bam.bai")
        tmp_summary = os.path.join(outdir, "web_summary.html")
        cell_count_job = cellranger.count(
            inputs=[readset.fastq1,readset.fastq2],
            output=tmp_bam,
            sample_id=readset.name + "_" + readset.sample_number,
            sample_name=readset.name,
            project=readset.project,
            ref_dir=self.get_reference_index()
        )

        
        output_bam = readset.bam + ".bam"
        output_bai = readset.bam + ".bai"
        output_summary = readset.bam + ".10x_summary.html"
        output_zip = readset.bam + "." + readset.run + "." + readset.lane_str + ".10x_outputs.zip"
        job = concat_jobs(
            [
                bash.mkdir(os.path.dirname(readset.bam + ".bam")),
                cell_count_job,
                bash.cp(
                    tmp_bam,
                    output_bam
                ),
                bash.cp(
                    tmp_bai,
                    output_bai
                ),
                bash.cp(
                    tmp_summary,
                    output_summary
                ),
                bash.zip(
                    outdir,
                    output_zip,
                    recursive=True
                )
            ],
            name="cellranger_count." + readset.name + "." + readset.run + "." + readset.lane,
            samples=[readset.sample],
            input_dependency=[readset.fastq1,readset.fastq2],
            output_dependency=[output_bam, output_bai, output_summary, output_zip]
        )
        return job

class VdjProcessingAligner(RNARunProcessingAligner):
    def get_annotation_files(self):
        return []

    def is_run_mark_duplicate(self):
        return False

    def get_reference_index(self):
        assembly_name = os.path.basename(self.genome_folder)
        ini_file = os.path.join(self.genome_folder, assembly_name + ".ini")
        if os.path.isfile(ini_file):
            config.parse_files([ini_file])

            if config.has_option("DEFAULT", "10x_vdj_transcriptome"):
                transcriptome = config.get("DEFAULT", "10x_vdj_transcriptome")
                trans_path = os.path.join(self.genome_folder, "genome", "10xGenomics", transcriptome)

                if trans_path is None or not os.path.isdir(trans_path):
                    """ Main transcriptome not found, we will try to resolve it if synomym exists """
                    if config.has_option("DEFAULT", "assembly_synonyms"):
                        synonym = config.get("DEFAULT", "assembly_synonyms")
                        trans_path = self.find_10x_synonym_reference(self.genome_folder, synonym, "10x_vdj_transcriptome")

            return trans_path
        else:
            return ""

    def get_reference_file(self):
        return os.path.join(self.get_reference_index(), "fasta", "genome.fa")

    def get_dictionary_file(self):
        folder_name = os.path.basename(self.genome_folder)
        return os.path.join(
            self.genome_folder,
            "genome",
            folder_name + ".dict"
        )

    def get_alignment_job(self, readset):
        outdir = os.path.join(
            self.output_dir,
            "10x_cellcount." + readset.lane,
            readset.name + "_" + readset.sample_number,
            "outs"
        )
        default_bam_name = "consensus"
        tmp_bam = os.path.join(outdir, readset.name + "." + default_bam_name + ".bam")
        tmp_bai = os.path.join(outdir, readset.name + "." + default_bam_name + ".bam.bai")
        tmp_summary = os.path.join(outdir, "web_summary.html")
        cell_vdj_job = cellranger.vdj(
            inputs=[readset.fastq1,readset.fastq2],
            output=tmp_bam,
            sample_id=readset.name + "_" + readset.sample_number,
            sample_name=readset.name,
            project=readset.project,
            ref_dir=self.get_reference_index()
        )

        output_bam = readset.bam + ".bam"
        output_bai = readset.bam + ".bai"
        output_summary = readset.bam + ".10x_summary.html"
        output_zip = readset.bam + "." + readset.run + "." + readset.lane + ".10x_outputs.zip"
        job = concat_jobs(
            [
                bash.mkdir(os.path.dirname(readset.bam + ".bam")),
                cell_vdj_job,
                bash.cp(
                    tmp_bam,
                    output_bam
                ),
                bash.cp(
                    tmp_bai,
                    output_bai
                ),
                bash.cp(
                    tmp_summary,
                    output_summary
                ),
                bash.zip(
                    outdir,
                    output_zip,
                    recursive=True
                )
            ],
            name="cellranger_vdj." + readset.name + "." + readset.run + "." + readset.lane,
            samples=[readset.sample],
            input_dependency=[readset.fastq1,readset.fastq2],
            output_dependency=[output_bam, output_bai, output_summary, output_zip]
        )
        return job

    def get_metrics_jobs(self, readset):
        jobs = []
        return jobs

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

    def get_dictionary_file(self):
        folder_name = os.path.basename(self.genome_folder)
        return os.path.join(
            self.genome_folder,
            "genome",
            folder_name + ".dict"
        )

    def get_alignment_job(self, readset):
        outdir = os.path.join(
            self.output_dir,
            "10x_cellcount." + readset.lane,
            readset.name + "_" + readset.sample_number,
            "outs"
        )
        default_bam_name = "possorted_bam"
        tmp_bam = os.path.join(outdir, readset.name + "." + default_bam_name + ".bam")
        tmp_bai = os.path.join(outdir, readset.name + "." + default_bam_name + ".bam.bai")
        tmp_summary = os.path.join(outdir, "web_summary.html")
        cell_atac_job = cellranger.atac(
            inputs=[readset.fastq1,readset.fastq2],
            output=tmp_bam,
            sample_id=readset.name + "_" + readset.sample_number,
            sample_name=readset.name,
            project=readset.project,
            ref_dir=self.get_reference_index()
        )

        output_bam = readset.bam + ".bam"
        output_bai = readset.bam + ".bai"
        output_summary = readset.bam + ".10x_summary.html"
        output_zip = readset.bam + "." + readset.run + "." + readset.lane + ".10x_outputs.zip"
        job = concat_jobs(
            [
                bash.mkdir(os.path.dirname(readset.bam + ".bam")),
                cell_atac_job,
                bash.cp(
                    tmp_bam,
                    output_bam
                ),
                bash.cp(
                    tmp_bai,
                    output_bai
                ),
                bash.cp(
                    tmp_summary,
                    output_summary
                ),
                bash.zip(
                    outdir,
                    output_zip,
                    recursive=True
                )
            ],
            name="cellranger_atac." + readset.name + "." + readset.run + "." + readset.lane,
            samples=[readset.sample],
            input_dependency=[readset.fastq1,readset.fastq2],
            output_dependency=[output_bam, output_bai, output_summary, output_zip]
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

