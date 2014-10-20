#!/usr/bin/env python

# Python Standard Modules
import logging
import math
import os
import re
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.readset import *
from bfx.sequence_dictionary import *

from bfx import bvatools
from bfx import bwa
from bfx import gatk
from bfx import gq_seq_utils
from bfx import igvtools
from bfx import metrics
from bfx import picard
from bfx import samtools
from bfx import snpeff
from bfx import tools
from bfx import vcftools
from pipelines import common

log = logging.getLogger(__name__)

class DnaSeq(common.Illumina):

    @property
    def sequence_dictionary(self):
        if not hasattr(self, "_sequence_dictionary"):
            self._sequence_dictionary = parse_sequence_dictionary_file(config.param('DEFAULT', 'genome_dictionary', type='filepath'))
        return self._sequence_dictionary

    def bwa_mem_picard_sort_sam(self):
        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_directory = os.path.join("alignment", readset.sample.name)
            readset_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam")

            if readset.run_type == "PAIRED_END":
                fastq1 = trim_file_prefix + "pair1.fastq.gz"
                fastq2 = trim_file_prefix + "pair2.fastq.gz"
            elif readset.run_type == "SINGLE_END":
                fastq1 = trim_file_prefix + "single.fastq.gz"
                fastq2 = None
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            job = concat_jobs([
                Job(command="mkdir -p " + os.path.dirname(readset_bam)),
                pipe_jobs([
                    bwa.mem(
                        fastq1,
                        fastq2,
                        read_group="'@RG" + \
                            "\tID:" + readset.name + \
                            "\tSM:" + readset.sample.name + \
                            ("\tLB:" + readset.library if readset.library else "") + \
                            ("\tPU:run" + readset.run + "_" + readset.lane if readset.run and readset.lane else "") + \
                            ("\tCN:" + config.param('bwa_mem', 'sequencing_center') if config.param('bwa_mem', 'sequencing_center', required=False) else "") + \
                            "\tPL:Illumina" + \
                            "'"
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        readset_bam,
                        "coordinate"
                    )
                ])
            ], name="bwa_mem_picard_sort_sam." + readset.name)

            # If this readset is unique for this sample, further BAM merging is not necessary.
            # Thus, create a sample BAM symlink to the readset BAM, along with its index.
            if len(readset.sample.readsets) == 1:
                readset_index = re.sub("\.bam$", ".bai", readset_bam)
                sample_bam = os.path.join(alignment_directory, readset.sample.name + ".sorted.bam")
                sample_index = re.sub("\.bam$", ".bai", sample_bam)
                job = concat_jobs([
                    job,
                    Job([readset_bam], [sample_bam], command="ln -s -f " + os.path.relpath(readset_bam, os.path.dirname(sample_bam)) + " " + sample_bam),
                    Job([readset_bam], [sample_index], command="ln -s -f " + os.path.relpath(readset_index, os.path.dirname(sample_index)) + " " + sample_index)
                ], name=job.name)

            jobs.append(job)
        return jobs

    def picard_merge_sam_files(self):
        jobs = []
        for sample in self.samples:
            # Skip samples with one readset only, since symlink has been created at align step
            if len(sample.readsets) > 1:
                alignment_directory = os.path.join("alignment", sample.name)
                inputs = [os.path.join(alignment_directory, readset.name + ".sorted.bam") for readset in sample.readsets]
                output = os.path.join(alignment_directory, sample.name + ".sorted.bam")

                job = picard.merge_sam_files(inputs, output)
                job.name = "picard_merge_sam_files." + sample.name
                jobs.append(job)
        return jobs

    def gatk_indel_realigner(self):
        jobs = []

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of realign jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            input = os.path.join(alignment_directory, sample.name + ".sorted.bam")

            if nb_jobs == 1:
                realign_prefix = os.path.join(realign_directory, "all")
                realign_intervals = realign_prefix + ".intervals"
                output_bam = realign_prefix + ".bam"
                sample_output_bam = os.path.join(alignment_directory, sample.name + ".realigned.qsorted.bam")
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + realign_directory),
                    gatk.realigner_target_creator(input, realign_intervals),
                    gatk.indel_realigner(input, output_bam, target_intervals=realign_intervals),
                    # Create sample realign symlink since no merging is required
                    Job([output_bam], [sample_output_bam], command="ln -s -f " + os.path.relpath(output_bam, os.path.dirname(sample_output_bam)) + " " + sample_output_bam)
                ], name="gatk_indel_realigner." + sample.name))

            else:
                # The first sequences are the longest to process.
                # Each of them must be processed in a separate job.
                unique_sequences_per_job = [sequence['name'] for sequence in self.sequence_dictionary[0:min(nb_jobs - 1, len(self.sequence_dictionary))]]

                # Create one separate job for each of the first sequences
                for sequence in unique_sequences_per_job:
                    realign_prefix = os.path.join(realign_directory, sequence)
                    realign_intervals = realign_prefix + ".intervals"
                    output_bam = realign_prefix + ".bam"
                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        Job(command="mkdir -p " + realign_directory),
                        gatk.realigner_target_creator(input, realign_intervals, intervals=[sequence]),
                        gatk.indel_realigner(input, output_bam, target_intervals=realign_intervals, intervals=[sequence])
                    ], name="gatk_indel_realigner." + sample.name + "." + sequence))

                # Create one last job to process the last remaining sequences and 'others' sequences
                realign_prefix = os.path.join(realign_directory, "others")
                realign_intervals = realign_prefix + ".intervals"
                output_bam = realign_prefix + ".bam"
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    Job(command="mkdir -p " + realign_directory),
                    gatk.realigner_target_creator(input, realign_intervals, exclude_intervals=unique_sequences_per_job),
                    gatk.indel_realigner(input, output_bam, target_intervals=realign_intervals, intervals=["unmapped"], exclude_intervals=unique_sequences_per_job)
                ], name="gatk_indel_realigner." + sample.name + ".others"))

        return jobs

    def merge_realigned(self):
        jobs = []

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', type='posint')

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            merged_realigned_bam = os.path.join(alignment_directory, sample.name + ".realigned.qsorted.bam")

            # if nb_jobs == 1, symlink has been created in indel_realigner and merging is not necessary
            if nb_jobs > 1:
                realigned_bams = [os.path.join(realign_directory, sequence['name'] + ".bam") for sequence in self.sequence_dictionary[0:min(nb_jobs - 1, len(self.sequence_dictionary))]]
                realigned_bams.append(os.path.join(realign_directory, "others.bam"))

                job = picard.merge_sam_files(realigned_bams, merged_realigned_bam)
                job.name = "merge_realigned." + sample.name
                jobs.append(job)

        return jobs

    def fix_mate_by_coordinate(self):
        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            input = alignment_file_prefix + "realigned.qsorted.bam"
            output_prefix = alignment_file_prefix + "matefixed.sorted"
            jobs.append(concat_jobs([
                bvatools.groupfixmate(input, output_prefix + ".tmp.bam"),
                samtools.sort(output_prefix + ".tmp.bam", output_prefix)
            ], name="fix_mate_by_coordinate." + sample.name))
        return jobs

    def picard_mark_duplicates(self):
        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            input = alignment_file_prefix + "matefixed.sorted.bam"
            output = alignment_file_prefix + "sorted.dup.bam"
            metrics_file = alignment_file_prefix + "sorted.dup.metrics"

            job = picard.mark_duplicates([input], output, metrics_file)
            job.name = "picard_mark_duplicates." + sample.name
            jobs.append(job)
        return jobs

    def recalibration(self):
        jobs = []
        for sample in self.samples:
            duplicate_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.")
            input = duplicate_file_prefix + "bam"
            print_reads_output = duplicate_file_prefix + "recal.bam"
            base_recalibrator_output = duplicate_file_prefix + "recalibration_report.grp"

            jobs.append(concat_jobs([
                gatk.base_recalibrator(input, base_recalibrator_output),
                gatk.print_reads(input, print_reads_output, base_recalibrator_output)
            ], name="recalibration." + sample.name))
        return jobs

    def metrics(self):
        jobs = []
        for sample in self.samples:
            recal_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.recal.")
            input = recal_file_prefix + "bam"

            job = picard.collect_multiple_metrics(input, recal_file_prefix + "all.metrics")
            job.name = "picard_collect_multiple_metrics." + sample.name
            jobs.append(job)

            # Compute genome coverage
            job = gatk.depth_of_coverage(input, recal_file_prefix + "all.coverage")
            job.name = "gatk_depth_of_coverage.genome." + sample.name
            jobs.append(job)

            # Compute CCDS coverage
            job = gatk.depth_of_coverage(input, recal_file_prefix + "CCDS.coverage", config.param('metrics', 'coverage_targets'))
            job.name = "gatk_depth_of_coverage.target." + sample.name
            jobs.append(job)

            job = bvatools.depth_of_coverage(
                input, 
                recal_file_prefix + "coverage.tsv", 
                bvatools.resolve_readset_coverage_bed(sample.readsets[0]), 
                config.param('bvatools_depth_of_coverage', 'other_options', required=False)
            )
            
            job.name = "bvatools_depth_of_coverage." + sample.name
            jobs.append(job)

            job = igvtools.compute_tdf(input, input + ".tdf")
            job.name = "igvtools_compute_tdf." + sample.name
            jobs.append(job)

            job = samtools.flagstat(input, recal_file_prefix + "bam.flagstat")
            job.name = "samtools_flagstat." + sample.name
            jobs.append(job)
        return jobs

    def picard_calculate_hs_metrics(self):
        jobs = []

        created_interval_lists = []

        for sample in self.samples:
            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            if coverage_bed:
                interval_list = re.sub("\.[^.]+$", ".interval_list", coverage_bed)

                recal_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.recal.")
                job = picard.calculate_hs_metrics(recal_file_prefix + "bam", recal_file_prefix + "onTarget.tsv", interval_list)
                if not interval_list in created_interval_lists:
                    job = concat_jobs([tools.bed2interval_list(None, coverage_bed, interval_list), job])
                    created_interval_lists.append(interval_list)

                job.name = "picard_calculate_hs_metrics." + sample.name
                jobs.append(job)
        return jobs

    def gatk_callable_loci(self):
        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")

            job = gatk.callable_loci(alignment_file_prefix + "sorted.dup.recal.bam", alignment_file_prefix + "callable.bed", alignment_file_prefix + "callable.summary.txt")
            job.name = "gatk_callable_loci." + sample.name
            jobs.append(job)

        return jobs

    def extract_common_snp_freq(self):
        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")

            job = bvatools.basefreq(alignment_file_prefix + "sorted.dup.recal.bam", alignment_file_prefix + "commonSNPs.alleleFreq.csv", config.param('extract_common_snp_freq', 'common_snp_positions', type='filepath'), 0)
            job.name = "extract_common_snp_freq." + sample.name
            jobs.append(job)

        return jobs

    def baf_plot(self):
        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")

            job = bvatools.ratiobaf(alignment_file_prefix + "commonSNPs.alleleFreq.csv", alignment_file_prefix + "ratioBAF", config.param('baf_plot', 'common_snp_positions', type='filepath'))
            job.name = "baf_plot." + sample.name
            jobs.append(job)

        return jobs

    def gatk_haplotype_caller(self):
        jobs = []

        nb_haplotype_jobs = config.param('gatk_haplotype_caller', 'nb_jobs', type='posint')
        if nb_haplotype_jobs > 50:
            log.warning("Number of haplotype jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            haplotype_directory = os.path.join(alignment_directory, "rawHaplotypeCaller")
            input = os.path.join(alignment_directory, sample.name + ".sorted.dup.recal.bam")

            if nb_haplotype_jobs == 1:
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    Job(command="mkdir -p " + haplotype_directory),
                    gatk.haplotype_caller(input, os.path.join(haplotype_directory, sample.name + ".hc.g.vcf"))
                ], name="gatk_haplotype_caller." + sample.name))

            else:
                unique_sequences_per_job = [sequence['name'] for sequence in self.sequence_dictionary[0:min(nb_haplotype_jobs - 1, len(self.sequence_dictionary))]]

                # Create one separate job for each of the first sequences
                for sequence in unique_sequences_per_job:
                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        Job(command="mkdir -p " + haplotype_directory),
                        gatk.haplotype_caller(input, os.path.join(haplotype_directory, sample.name + "." + sequence + ".hc.g.vcf"), intervals=[sequence])
                    ], name="gatk_haplotype_caller." + sample.name + "." + sequence))

                # Create one last job to process the last remaining sequences and 'others' sequences
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    Job(command="mkdir -p " + haplotype_directory),
                    gatk.haplotype_caller(input, os.path.join(haplotype_directory, sample.name + ".others.hc.g.vcf"), exclude_intervals=unique_sequences_per_job)
                ], name="gatk_haplotype_caller." + sample.name + ".others"))

        return jobs

    def merge_and_call_gvcf(self):
        jobs = []
        nb_haplotype_jobs = config.param('gatk_haplotype_caller', 'nb_jobs', type='posint')

        for sample in self.samples:
            haplotype_file_prefix = os.path.join("alignment", sample.name, "rawHaplotypeCaller", sample.name)
            output_haplotype_file_prefix = os.path.join("alignment", sample.name, sample.name)
            if nb_haplotype_jobs == 1:
                gvcfs_to_merge = [haplotype_file_prefix + ".hc.g.vcf"]
            else:
                gvcfs_to_merge = [haplotype_file_prefix + "." + sequence['name'] + ".hc.g.vcf" for sequence in self.sequence_dictionary[0:min(nb_haplotype_jobs - 1, len(self.sequence_dictionary))]]
                gvcfs_to_merge.append(haplotype_file_prefix + ".others.hc.g.vcf")

            jobs.append(concat_jobs([
                gatk.cat_variants(gvcfs_to_merge, output_haplotype_file_prefix + ".hc.g.vcf"),
                gatk.genotype_gvcfs([output_haplotype_file_prefix + ".hc.g.vcf"], output_haplotype_file_prefix + ".hc.vcf")
            ], name="merge_and_call_gvcf." + sample.name))

        return jobs

    def dna_sample_metrics(self):
        job = concat_jobs([
            Job(command="mkdir -p metrics"),
            metrics.dna_sample_metrics("alignment", "metrics/SampleMetrics.stats", config.param('DEFAULT', 'experiment_type'))
        ], name="dna_sample_metrics")
        job.input_files = [os.path.join("alignment", sample.name, sample.name + ".sorted.dup.metrics") for sample in self.samples]
        return [job]

    def generate_approximate_windows(self, nb_jobs):
        if nb_jobs <= len(self.sequence_dictionary):
            return [sequence['name'] + ":1-" + str(sequence['length']) for sequence in self.sequence_dictionary]
        else:
            total_length = sum([sequence['length'] for sequence in self.sequence_dictionary])
            approximate_window_size = int(math.floor(total_length / (nb_jobs - len(self.sequence_dictionary))))
            windows = []

            for sequence in self.sequence_dictionary:
                for start, end in [[pos, min(pos + approximate_window_size - 1, sequence['length'])] for pos in range(1, sequence['length'] + 1, approximate_window_size)]:
                    windows.append(sequence['name'] + ":" + str(start) + "-" + str(end))
            return windows

    def rawmpileup(self):
        jobs = []
        for sample in self.samples:
            mpileup_directory = os.path.join("alignment", sample.name, "mpileup")

            for sequence in self.sequence_dictionary:
                output = os.path.join(mpileup_directory, sample.name + "." + sequence['name'] + ".mpileup.gz")
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + mpileup_directory),
                    pipe_jobs([
                        samtools.mpileup([os.path.join("alignment", sample.name, sample.name + ".sorted.dup.recal.bam")], None, config.param('rawmpileup', 'mpileup_other_options'), sequence['name']),
                        Job(output_files=[output], command="gzip -1 -c > " + output)
                    ])], name="rawmpileup." + sample.name + "." + sequence['name']))

        return jobs

    def rawmpileup_cat(self):
        jobs = []
        for sample in self.samples:
            mpileup_file_prefix = os.path.join("alignment", sample.name, "mpileup", sample.name + ".")
            mpileup_inputs = [mpileup_file_prefix + sequence['name'] + ".mpileup.gz" for sequence in self.sequence_dictionary]

            gzip_output = mpileup_file_prefix + "mpileup.gz"
            job = Job(mpileup_inputs, [gzip_output])
            job.command = "zcat \\\n  " + " \\\n  ".join(mpileup_inputs) + " | \\\n  gzip -c --best > " + gzip_output
            job.name = "rawmpileup_cat." + sample.name
            jobs.append(job)
        return jobs

    def snp_and_indel_bcf(self):
        jobs = []
        input_bams = [os.path.join("alignment", sample.name, sample.name + ".sorted.dup.recal.bam") for sample in self.samples]
        nb_jobs = config.param('snp_and_indel_bcf', 'approximate_nb_jobs', type='posint')
        output_directory = "variants/rawBCF"
        bcftools_view_options = "-bvcg"

        if nb_jobs == 1:
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + output_directory),
                pipe_jobs([
                    samtools.mpileup(input_bams, None, config.param('snp_and_indel_bcf', 'mpileup_other_options')),
                    samtools.bcftools_view("-", os.path.join(output_directory, "allSamples.bcf"), bcftools_view_options),
                ])], name="snp_and_indel_bcf.allSamples"))

        else:
            for region in self.generate_approximate_windows(nb_jobs):
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + output_directory),
                    pipe_jobs([
                        samtools.mpileup(input_bams, None, config.param('snp_and_indel_bcf', 'mpileup_other_options'), region),
                        samtools.bcftools_view("-", os.path.join(output_directory, "allSamples." + region + ".bcf"), bcftools_view_options),
                    ])], name="snp_and_indel_bcf.allSamples." + re.sub(":", "_", region)))

        return jobs

    def merge_filter_bcf(self):
        nb_jobs = config.param('snp_and_indel_bcf', 'approximate_nb_jobs', type='posint')

        if nb_jobs == 1:
            inputs = ["variants/rawBCF/allSamples.bcf"]
        else:
            inputs = ["variants/rawBCF/allSamples." + region + ".bcf" for region in self.generate_approximate_windows(nb_jobs)]
        output_file_prefix = "variants/allSamples.merged."

        bcf = output_file_prefix + "bcf"
        job = concat_jobs([
            samtools.bcftools_cat(inputs, bcf),
            samtools.bcftools_view(bcf, output_file_prefix + "flt.vcf")
        ])
        job.name = "merge_filter_bcf"
        return [job]

    def filter_nstretches(self):
        job = tools.filter_long_indel("variants/allSamples.merged.flt.vcf", "variants/allSamples.merged.flt.NFiltered.vcf")
        job.name = "filter_nstretches"
        return [job]

    def flag_mappability(self):
        job = vcftools.annotate_mappability("variants/allSamples.merged.flt.NFiltered.vcf", "variants/allSamples.merged.flt.mil.vcf")
        job.name = "flag_mappability"
        return [job]

    def snp_id_annotation(self):
        job = snpeff.snpsift_annotate("variants/allSamples.merged.flt.mil.vcf", "variants/allSamples.merged.flt.mil.snpId.vcf")
        job.name = "snp_id_annotation"
        return [job]

    def snp_effect(self):
        job = snpeff.compute_effects("variants/allSamples.merged.flt.mil.snpId.vcf", "variants/allSamples.merged.flt.mil.snpId.snpeff.vcf", split=True)
        job.name = "snp_effect"
        return [job]

    def dbnsfp_annotation(self):
        job = snpeff.snpsift_dbnsfp("variants/allSamples.merged.flt.mil.snpId.snpeff.vcf", "variants/allSamples.merged.flt.mil.snpId.snpeff.dbnsfp.vcf")
        job.name = "dbnsfp_annotation"
        return [job]

    def metrics_vcf_stats(self):
        variants_file_prefix = "variants/allSamples.merged.flt.mil.snpId."

        job = metrics.vcf_stats(variants_file_prefix + "vcf", variants_file_prefix + "snpeff.vcf.part_changeRate.tsv", variants_file_prefix + "snpeff.vcf.statsFile.txt")
        job.name = "metrics_change_rate"
        return [job]


    def metrics_snv_graph_metrics(self):
        variants_file_prefix = "variants/allSamples.merged.flt.mil.snpId."
        job = metrics.snv_graph_metrics(variants_file_prefix + "snpeff.vcf.statsFile.txt", "metrics/allSamples.SNV")
        job.output_files = ["metrics/allSamples.SNV.SummaryTable.tsv"]
        job.name = "metrics_snv_graph"

        return [job]

    def deliverable(self):
        job = gq_seq_utils.report(os.path.abspath(config.filepath), self.output_dir, "DNAseq", "deliverables")
        job.input_files = [
            "metrics/SampleMetrics.stats",
            "variants/allSamples.merged.flt.vcf",
            "metrics/allSamples.SNV.SummaryTable.tsv"
        ]
        job.name = "deliverable"
        return [job]

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.bwa_mem_picard_sort_sam,
            self.picard_merge_sam_files,
            self.gatk_indel_realigner,
            self.merge_realigned,
            self.fix_mate_by_coordinate,
            self.picard_mark_duplicates,
            self.recalibration,
            self.metrics,
            self.picard_calculate_hs_metrics,
            self.gatk_callable_loci,
            self.extract_common_snp_freq,
            self.baf_plot,
            self.gatk_haplotype_caller,
            self.merge_and_call_gvcf,
            self.dna_sample_metrics,
            self.rawmpileup,
            self.rawmpileup_cat,
            self.snp_and_indel_bcf,
            self.merge_filter_bcf,
            self.filter_nstretches,
            self.flag_mappability,
            self.snp_id_annotation,
            self.snp_effect,
            self.dbnsfp_annotation,
            self.metrics_vcf_stats,
            self.metrics_snv_graph_metrics,
            self.deliverable
        ]

if __name__ == '__main__': 
    DnaSeq().submit_jobs()
