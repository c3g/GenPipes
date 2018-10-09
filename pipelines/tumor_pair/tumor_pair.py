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

# Python Standard Modules
import logging
import math
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.sample_tumor_pairs import *
from bfx.sequence_dictionary import *

from bfx import picard
from bfx import sambamba
from bfx import samtools
from bfx import vcflib
from bfx import htslib
from bfx import varscan
from bfx import bamreadcount
from bfx import bcftools
from bfx import gq_seq_utils
from bfx import gatk
from bfx import tools
from bfx import bed_file
from bfx import vardict
from bfx import bcbio_variation_recall
from bfx import vt
from bfx import snpeff
from bfx import gemini
from bfx import conpair
from pipelines.dnaseq import dnaseq

log = logging.getLogger(__name__)

class TumorPair(dnaseq.DnaSeq):
    """
    Tumor Pair Pipeline
    =================

    The Tumor Pair pipeline inherits the initial bam preparation steps of the DNA-Seq pipeline with the exception of the
    indel realignment (IR) step. In the tumor pipeline the IR step utilizes both the normal and tumor bam to further reduce
    false positives (FPs) in and around indels. The tumor pipeline deviates from the DNA-seq pipeline at the variant calling step.
    At this point, a paired caller is used to call SNVs and Indels from the pairs given as input. Additional, muliple cancer callers
    are utilized using an ensemble approach and SNVs and Indels seen in at least 2 different callers are retained for further
    investigation.

    Example command:
    python tumor_pair.py -c a.ini b.base.ini -s x-y,z -r readset.tsv -p pairs.csv

    -c ini files: multiple can be specified e.g WGS or exome, or different clusters e.g. base (abacus) or guillimin

    -r readset: derived from GQ lims or made yourself. See : https://bitbucket.org/mugqic/mugqic_pipelines#markdown-header-readset-file

    -p pairs : format - patient_name,normal_sample_name,tumor_sample_name
    """

    def __init__(self, protocol=None):
        self._protocol=protocol
        # Add pipeline specific arguments
        self.argparser.add_argument("-p", "--pairs", help="pairs file", type=file)
        super(TumorPair, self).__init__(protocol)

    @property
    def tumor_pairs(self):
        if not hasattr(self, "_tumor_pairs"):
            self._tumor_pairs = parse_tumor_pair_file(self.args.pairs.name, self.samples)
        return self._tumor_pairs

    def sequence_dictionary_variant(self):
        if not hasattr(self, "_sequence_dictionary_variant"):
            self._sequence_dictionary_variant = parse_sequence_dictionary_file(config.param('DEFAULT', 'genome_dictionary', type='filepath'), variant=True)
        return self._sequence_dictionary_variant

    def generate_approximate_windows(self, nb_jobs):
        if nb_jobs <= len(self.sequence_dictionary_variant()):
            return [sequence['name'] + ":1-" + str(sequence['length']) for sequence in self.sequence_dictionary_variant()]
        else:
            total_length = sum([sequence['length'] for sequence in self.sequence_dictionary_variant()])
            approximate_window_size = int(math.floor(total_length / (nb_jobs - len(self.sequence_dictionary_variant()))))
            windows = []

            for sequence in self.sequence_dictionary_variant():
                for start, end in [[pos, min(pos + approximate_window_size - 1, sequence['length'])] for pos in range(1, sequence['length'] + 1, approximate_window_size)]:
                    windows.append(sequence['name'] + ":" + str(start) + "-" + str(end))
            return windows

    def sambamba_merge_sam_files(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [Picard](http://broadinstitute.github.io/picard/).

        This step takes as input files:

        1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
        2. Else, BAM files from the readset file
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            # Find input readset BAMs first from previous bwa_mem_picard_sort_sam job, then from original BAMs in the readset sheet.
            readset_bams = self.select_input_files([[os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam") for readset in sample.readsets], [readset.bam for readset in sample.readsets]])
            sample_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")

            mkdir_job = Job(command="mkdir -p " + os.path.dirname(sample_bam), samples=[sample])

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            if len(sample.readsets) == 1:
                readset_bam = readset_bams[0]
                if os.path.isabs(readset_bam):
                    target_readset_bam = readset_bam
                else:
                    target_readset_bam = os.path.relpath(readset_bam, alignment_directory)
                readset_index = re.sub("\.bam$", ".bai", readset_bam)
                target_readset_index = re.sub("\.bam$", ".bai", target_readset_bam)
                sample_index = re.sub("\.bam$", ".bai", sample_bam)

                job = concat_jobs([
                    mkdir_job,
                    Job([readset_bam], [sample_bam], command="ln -s -f " + target_readset_bam + " " + sample_bam, removable_files=[sample_bam]),
                    Job([readset_index], [sample_index], command="ln -s -f " + target_readset_index + " " + sample_index, removable_files=[sample_index])
                ], name="symlink_readset_sample_bam." + sample.name)

            elif len(sample.readsets) > 1:
                job = concat_jobs([
                    mkdir_job,
                    sambamba.merge(readset_bams, sample_bam)
                ])
                job.name = "sambamba_merge_sam_files." + sample.name

            jobs.append(job)

        return jobs


    def gatk_indel_realigner(self):
        """
        Insertion and deletion realignment is performed on regions where multiple base mismatches
        are preferred over indels by the aligner since it can appear to be less costly by the algorithm.
        Such regions will introduce false positive variant calls which may be filtered out by realigning
        those regions properly. Realignment is done using [GATK](https://www.broadinstitute.org/gatk/).
        The reference genome is divided by a number regions given by the `nb_jobs` parameter.

        Note: modified to use both normal and tumor bams to reduce FPs around indels

        """

        jobs = []

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of realign jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            normal_alignment_directory = os.path.join("alignment", tumor_pair.normal.name)
            normal_realign_directory = os.path.join(normal_alignment_directory, "realign")
            tumor_alignment_directory = os.path.join("alignment", tumor_pair.tumor.name)
            tumor_realign_directory = os.path.join(tumor_alignment_directory, "realign")

            input_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")
            input_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")

            if nb_jobs == 1:
                realign_intervals = os.path.join(tumor_realign_directory, "all.intervals")
                bam_postfix = ".all.realigned.bam"
                normal_bam = os.path.join(tumor_pair.normal.name + ".sorted.all.realigned.bam")
                normal_index = re.sub("\.bam$", ".bai", normal_bam)
                normal_output_bam = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".realigned.qsorted.bam")
                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                tumor_bam = os.path.join( tumor_pair.tumor.name + ".sorted.all.realigned.bam")
                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                tumor_output_bam = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".realigned.qsorted.bam")
                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + normal_realign_directory, removable_files=[normal_realign_directory], samples=[tumor_pair.normal]),
                    Job(command="mkdir -p " + tumor_realign_directory, removable_files=[tumor_realign_directory], samples=[tumor_pair.tumor]),
                    gatk.realigner_target_creator(input_normal, realign_intervals, input2=input_tumor),
                    gatk.indel_realigner(input_normal, input2=input_tumor, output_norm_dep=normal_output_bam, output_tum_dep=tumor_output_bam, target_intervals=realign_intervals, optional=bam_postfix),
                    # Move sample realign
                    Job([input_normal], [normal_output_bam], command="mv " + normal_bam + " " + normal_output_bam + " && mv " + normal_index + " " + normal_output_index),
                    Job([input_tumor], [tumor_output_bam], command="mv " + tumor_bam + " " + tumor_output_bam + " && mv " + tumor_index + " " + tumor_output_index)
                ], name="gatk_indel_realigner." + tumor_pair.name))

            else:
                # The first sequences are the longest to process.
                # Each of them must be processed in a separate job.
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                # Create one separate job for each of the first sequences
                for idx,sequences in enumerate(unique_sequences_per_job):
                    realign_prefix = os.path.join(tumor_realign_directory, str(idx))
                    realign_intervals = realign_prefix + ".intervals"
                    intervals = sequences
                    if str(idx) == 0:
                        intervals.append("unmapped")
                    bam_postfix = ".realigned." + str(idx) + ".bam"
                    normal_bam = os.path.join(tumor_pair.normal.name + ".sorted.realigned." + str(idx)  + ".bam")
                    normal_index = re.sub("\.bam$", ".bai", normal_bam)
                    tumor_bam = os.path.join(tumor_pair.tumor.name + ".sorted.realigned." + str(idx)  + ".bam")
                    tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                    normal_output_bam = os.path.join(normal_realign_directory, tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam")
                    normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                    tumor_output_bam = os.path.join(tumor_realign_directory, tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam")
                    tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        Job(command="mkdir -p " + normal_realign_directory, removable_files=[normal_realign_directory], samples=[tumor_pair.normal]),
                        Job(command="mkdir -p " + tumor_realign_directory, removable_files=[tumor_realign_directory], samples=[tumor_pair.tumor]),
                        gatk.realigner_target_creator(input_normal, realign_intervals, input2=input_tumor, intervals=intervals),
                        gatk.indel_realigner(input_normal, input2=input_tumor, output_norm_dep=normal_output_bam, output_tum_dep=tumor_output_bam, target_intervals=realign_intervals, intervals=intervals, optional=bam_postfix),
                        Job([input_normal], [normal_output_bam], command="mv " + normal_bam + " " + normal_output_bam + " && mv " + normal_index + " " + normal_output_index),
                        Job([input_tumor], [tumor_output_bam], command="mv " + tumor_bam + " " + tumor_output_bam + " && mv " + tumor_index + " " + tumor_output_index)
                    ], name="gatk_indel_realigner." + tumor_pair.name + "." + str(idx)))

                # Create one last job to process the last remaining sequences and 'others' sequences
                realign_prefix = os.path.join(tumor_realign_directory, "others")
                realign_intervals = realign_prefix + ".intervals"
                bam_postfix = ".realigned.others.bam"
                normal_bam = os.path.join(tumor_pair.normal.name + ".sorted.realigned.others.bam")
                normal_index = re.sub("\.bam$", ".bai", normal_bam)
                tumor_bam = os.path.join(tumor_pair.tumor.name + ".sorted.realigned.others.bam")
                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                normal_output_bam = os.path.join(normal_realign_directory, tumor_pair.normal.name + ".sorted.realigned.others.bam")
                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                tumor_output_bam = os.path.join(tumor_realign_directory, tumor_pair.tumor.name + ".sorted.realigned.others.bam")
                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    Job(command="mkdir -p " + normal_realign_directory, removable_files=[normal_realign_directory], samples=[tumor_pair.normal]),
                    Job(command="mkdir -p " + tumor_realign_directory, removable_files=[tumor_realign_directory], samples=[tumor_pair.tumor]),
                    gatk.realigner_target_creator(input_normal, realign_intervals, input2=input_tumor, exclude_intervals=unique_sequences_per_job_others),
                    gatk.indel_realigner(input_normal, input2=input_tumor, output_norm_dep=normal_output_bam, output_tum_dep=tumor_output_bam, target_intervals=realign_intervals, exclude_intervals=unique_sequences_per_job_others, optional=bam_postfix),
                    Job([input_normal], [normal_output_bam], command="mv " + normal_bam + " " + normal_output_bam + " && mv " + normal_index + " " + normal_output_index),
                    Job([input_tumor], [tumor_output_bam], command="mv " + tumor_bam + " " + tumor_output_bam + " && mv " + tumor_index + " " + tumor_output_index)
                ], name="gatk_indel_realigner." + tumor_pair.name + ".others"))

        return jobs

    def merge_realigned(self):
        """
        BAM files of regions of realigned reads are merged per sample using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', type='posint')

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            merged_realigned_bam = os.path.join(alignment_directory, sample.name + ".realigned.qsorted.bam")

            if nb_jobs > 1:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                inputBAMs = []
                for idx,sequences in enumerate(unique_sequences_per_job):
                    inputBAMs.append(os.path.join(realign_directory, sample.name + ".sorted.realigned." + str(idx) + ".bam"))
                inputBAMs.append(os.path.join(realign_directory, sample.name + ".sorted.realigned.others.bam"))

                job = picard.merge_sam_files(inputBAMs, merged_realigned_bam)
                job.name = "merge_realigned." + sample.name
                job.samples = [sample]
                jobs.append(job)

        report_file = os.path.join("report", "DnaSeq.gatk_indel_realigner.md")
        jobs.append(
            Job(
                [os.path.join("alignment", sample.name, sample.name + ".realigned.sorted.bam") for sample in self.samples],
                [report_file],
                command="""\
mkdir -p report && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="merge_realigned_report")
        )

        return jobs

    def sambamba_merge_realigned(self):
        """
        BAM files of regions of realigned reads are merged per sample using [Sambamba](http://lomereiter.github.io/sambamba/index.html).
        """

        jobs = []

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', type='posint')

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            merged_realigned_bam = os.path.join(alignment_directory, sample.name + ".realigned.sorted.bam")

            # if nb_jobs == 1, symlink has been created in indel_realigner and merging is not necessary
            if nb_jobs > 1:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                inputBAMs = []
                for idx,sequences in enumerate(unique_sequences_per_job):
                    inputBAMs.append(os.path.join(realign_directory, sample.name + ".sorted.realigned." + str(idx) + ".bam"))
                inputBAMs.append(os.path.join(realign_directory, sample.name + ".sorted.realigned.others.bam"))
                #realigned_bams = [os.path.join(realign_directory, sequence['name'] + ".bam") for sequence in self.sequence_dictionary[0:min(nb_jobs - 1, len(self.sequence_dictionary))]]
                #realigned_bams.append(os.path.join(realign_directory, "others.bam"))

                job = sambamba.merge(inputBAMs, merged_realigned_bam)
                job.name = "sambamba_merge_realigned." + sample.name
                job.samples = [sample]
                jobs.append(job)

        report_file = os.path.join("report", "DnaSeq.gatk_indel_realigner.md")
        jobs.append(
            Job(
                [os.path.join("alignment", sample.name, sample.name + ".realigned.sorted.bam") for sample in self.samples],
                [report_file],
                command="""\
mkdir -p report && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="merge_realigned_report")
        )

        return jobs

    def sambamba_mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            input = alignment_file_prefix + "realigned.sorted.bam"
            output = alignment_file_prefix + "sorted.dup.bam"
            metrics_file = alignment_file_prefix + "sorted.dup.metrics"

            job = sambamba.markdup(input, output, os.path.join("alignment", sample.name, sample.name))
            job.name = "sambamba_mark_duplicates." + sample.name
            job.samples = [sample]
            jobs.append(job)

        report_file = os.path.join("report", "DnaSeq.picard_mark_duplicates.md")
        jobs.append(
            Job(
                [os.path.join("alignment", sample.name, sample.name + ".sorted.dup.bam") for sample in self.samples],
                [report_file],
                command="""\
mkdir -p report && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="picard_mark_duplicates_report")
        )

        return jobs

    def conpair_concordance_contamination(self):
        """

        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            metrics_directory = os.path.join("metrics")
            input_normal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam")
            input_tumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")
            pileup_normal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".gatkPileup")
            pileup_tumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".gatkPileup")

            concordance_out = os.path.join(metrics_directory, tumor_pair.name + ".concordance.tsv")
            contamination_out = os.path.join(metrics_directory, tumor_pair.name + ".contamination.tsv")

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal]),
                conpair.pileup(input_normal, pileup_normal),
            ], name="conpair_concordance_contamination.pileup." + tumor_pair.normal.name))

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.tumor]),
                conpair.pileup(input_tumor, pileup_tumor),
            ], name="conpair_concordance_contamination.pileup." + tumor_pair.tumor.name))

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + metrics_directory, samples=[tumor_pair.normal, tumor_pair.tumor]),
                conpair.concordance(pileup_normal, pileup_tumor, concordance_out),
                conpair.contamination(pileup_normal, pileup_tumor, contamination_out)
            ], name="conpair_concordance_contamination." + tumor_pair.name))

        return jobs


    def rawmpileup_panel(self):
        """
        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
        One packaged mpileup file is created per sample/chromosome.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            bedfile=config.param('rawmpileup_panel', 'panel')

            for sequence in self.sequence_dictionary_variant():
                normal_output = os.path.join(varscan_directory, tumor_pair.normal.name + "." + sequence['name'] + ".mpileup")
                tumor_output = os.path.join(varscan_directory, tumor_pair.tumor.name + "." + sequence['name'] + ".mpileup")
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory], samples=[tumor_pair.normal, tumor_pair.tumor]),
                    samtools.mpileup([os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam")], normal_output, config.param('rawmpileup_panel', 'mpileup_other_options'), region=sequence['name'], regionFile=bedfile),
                    samtools.mpileup([os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")], tumor_output, config.param('rawmpileup_panel', 'mpileup_other_options'), region=sequence['name'], regionFile=bedfile),
                ], name="rawmpileup_panel." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def paired_varscan2_panel(self):
        """
        Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name , "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            for sequence in self.sequence_dictionary_variant():
                input_normal = os.path.join(varscan_directory, tumor_pair.normal.name + "." + sequence['name'] + ".mpileup")
                input_tumor = os.path.join(varscan_directory, tumor_pair.tumor.name + "." + sequence['name'] + ".mpileup")

                output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'])
                output_snp =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf")
                output_indel =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf")
                output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + ".varscan2." + sequence['name'] + ".vcf.gz")

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory], samples=[tumor_pair.normal, tumor_pair.tumor]),
                    varscan.somatic(input_normal, input_tumor, output, config.param('varscan2_somatic_panel', 'other_options'), output_vcf_dep=output_vcf_gz, output_snp_dep=output_snp, output_indel_dep=output_indel),
                    htslib.bgzip_tabix_vcf(output_snp, os.path.join(varscan_directory, tumor_pair.name + ".snp." + sequence['name'] + ".vcf.gz")),
                    htslib.bgzip_tabix_vcf(output_indel, os.path.join(varscan_directory, tumor_pair.name + ".indel." + sequence['name'] + ".vcf.gz")),
                    pipe_jobs([
                        bcftools.concat([os.path.join(varscan_directory, tumor_pair.name + ".snp." + sequence['name'] + ".vcf.gz"), os.path.join(varscan_directory, tumor_pair.name + ".indel." + sequence['name'] + ".vcf.gz")], None),
                        Job([None], [None], command="sed 's/TUMOR/"+ tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"+ tumor_pair.normal.name + "/g' " ),
                        htslib.bgzip_tabix_vcf(None, output_vcf_gz),
                    ]),
                ], name = "varscan2_somatic_panel." + tumor_pair.name + "." + sequence['name'] ) )

        return jobs

    def merge_varscan2_panel(self):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            all_inputs = [os.path.join(varscan_directory, tumor_pair.name + ".varscan2." + sequence['name'] + ".vcf.gz") for sequence in self.sequence_dictionary_variant()]

            all_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz")
            somatic_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz")
            germline_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vcf.gz")

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    bcftools.concat(all_inputs, None),
                    tools.fix_varscan_output(None, None),
                    Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                    Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                    Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),
                    htslib.bgzip_tabix_vcf(None, all_output),
                ]),
                pipe_jobs([
                    bcftools.view(all_output, None, config.param('merge_varscan2', 'somatic_filter_options')),
                    htslib.bgzip_tabix_vcf(None, somatic_output),
                ]),
                pipe_jobs([
                    bcftools.view(all_output, None, config.param('merge_varscan2', 'germline_loh_filter_options')),
                    htslib.bgzip_tabix_vcf(None, germline_output),
                ]),
            ], name = "merge_varscan2." + tumor_pair.name ))

        return jobs

    def preprocess_vcf_panel(self):
        """
        Preprocess vcf for loading into a annotation database - gemini : http://gemini.readthedocs.org/en/latest/index.html
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt) and
        vcf FORMAT modification for correct loading into gemini
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            prefix = os.path.join(pair_directory, tumor_pair.name)
            outputPreprocess = prefix + ".prep.vt.vcf.gz"
            outputFix = prefix + ".varscan2.somatic.vt.vcf.gz"

            outputPrepGermline = prefix + ".germline_loh.prep.vt.vcf.gz"
            outputFixGermline = prefix + ".varscan2.germline_loh.vt.vcf.gz"

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    vt.decompose_and_normalize_mnps( prefix + ".varscan2.somatic.vcf.gz" , None),
                    htslib.bgzip_tabix_vcf(None, prefix + ".prep.vt.vcf.gz"),
                ]),
                tools.preprocess_varscan( prefix + ".prep.vt.vcf.gz",  outputFix ),
                #Job([outputPreprocess], [outputFix], command="zcat " + outputPreprocess + " | grep -v 'ID=AD_O' | awk ' BEGIN {OFS=\"\\t\"; FS=\"\\t\"} {if (NF > 8) {for (i=9;i<=NF;i++) {x=split($i,na,\":\") ; if (x > 1) {tmp=na[1] ; for (j=2;j<x;j++){if (na[j] == \"AD_O\") {na[j]=\"AD\"} ; if (na[j] != \".\") {tmp=tmp\":\"na[j]}};$i=tmp}}};print $0} ' | bgzip -cf >  " + outputFix),
            ], name="preprocess_vcf_panel.somatic." + tumor_pair.name ))

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    vt.decompose_and_normalize_mnps( prefix + ".varscan2.germline_loh.vcf.gz" , None),
                    htslib.bgzip_tabix_vcf(None,  prefix + ".germline_loh.prep.vt.vcf.gz"),
                ]),
                tools.preprocess_varscan( prefix + ".germline_loh.prep.vt.vcf.gz",  outputFixGermline ),
                #Job([outputPrepGermline], [outputFixGermline], command="zcat " + outputPrepGermline + " | grep -v 'ID=AD_O' | awk ' BEGIN {OFS=\"\\t\"; FS=\"\\t\"} {if (NF > 8) {for (i=9;i<=NF;i++) {x=split($i,na,\":\") ; if (x > 1) {tmp=na[1] ; for (j=2;j<x;j++){if (na[j] == \"AD_O\") {na[j]=\"AD\"} ; if (na[j] != \".\") {tmp=tmp\":\"na[j]}};$i=tmp}}};print $0} ' | bgzip -cf >  " + outputFixGermline),
            ], name="preprocess_vcf_panel.germline." + tumor_pair.name ))

        return jobs

    def snp_effect_panel(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            if not os.path.exists(varscan_directory):
                os.makedirs(varscan_directory)

            input_somatic = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            output_somatic = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.snpeff.vcf")
            output_somatic_gz = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.snpeff.vcf.gz")

            input_germline = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vt.vcf.gz")
            output_germline = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vt.snpeff.vcf")
            output_germline_gz = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vt.snpeff.vcf.gz")

            cancer_pair_filename = os.path.join(varscan_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write( tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n" )

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                snpeff.compute_effects(input_somatic, output_somatic, cancer_sample_file=cancer_pair_filename, options=config.param('compute_cancer_effects_somatic', 'options')),
                htslib.bgzip_tabix_vcf(output_somatic, output_somatic_gz),
            ], name = "compute_cancer_effects_somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                snpeff.compute_effects(input_germline, output_germline, cancer_sample_file=cancer_pair_filename, options=config.param('compute_cancer_effects_germline', 'options')),
                htslib.bgzip_tabix_vcf(output_germline, output_germline_gz),
            ], name = "compute_cancer_effects_germline." + tumor_pair.name))

        return jobs

    def gemini_annotations_panel(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        gemini_module=config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2],gemini_module[-1]])

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            if not os.path.exists(varscan_directory):
                os.makedirs(varscan_directory)

            temp_dir = os.path.join(os.getcwd(), pair_directory)
            gemini_prefix = os.path.join(pair_directory, tumor_pair.name)

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                gemini.gemini_annotations( gemini_prefix + ".varscan2.somatic.vt.snpeff.vcf.gz", gemini_prefix + ".somatic.gemini." + gemini_version + ".db", temp_dir)
            ], name="gemini_annotations.somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                gemini.gemini_annotations( gemini_prefix + ".varscan2.germline_loh.vt.snpeff.vcf.gz", gemini_prefix + ".germline_loh.gemini." + gemini_version + ".db", temp_dir)
            ], name="gemini_annotations.germline." + tumor_pair.name))

        return jobs


    def rawmpileup(self):
        """
        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
        One packaged mpileup file is created per sample/chromosome.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            for sequence in self.sequence_dictionary_variant():
                normal_output = os.path.join(varscan_directory, tumor_pair.normal.name + "." + sequence['name'] + ".mpileup")
                tumor_output = os.path.join(varscan_directory, tumor_pair.tumor.name + "." + sequence['name'] + ".mpileup")
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory], samples=[tumor_pair.normal, tumor_pair.tumor]),
                    samtools.mpileup([os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")], normal_output, config.param('rawmpileup', 'mpileup_other_options'), sequence['name']),
                    samtools.mpileup([os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")], tumor_output, config.param('rawmpileup', 'mpileup_other_options'), sequence['name']),
                ], name="rawmpileup." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def rawmpileup_cat(self):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            mpileup_normal_file_prefix = os.path.join(varscan_directory, tumor_pair.normal.name + ".")
            mpileup_normal_inputs = [mpileup_normal_file_prefix + sequence['name'] + ".mpileup" for sequence in self.sequence_dictionary_variant()]

            mpileup_tumor_file_prefix = os.path.join(varscan_directory, tumor_pair.tumor.name + ".")
            mpileup_tumor_inputs = [mpileup_tumor_file_prefix + sequence['name'] + ".mpileup" for sequence in self.sequence_dictionary_variant()]

            normal_output = mpileup_normal_file_prefix + "mpileup"
            tumor_output = mpileup_tumor_file_prefix + "mpileup"

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory]),
                Job(mpileup_normal_inputs, [normal_output], command = "cat \\\n  " + " \\\n  ".join(mpileup_normal_inputs) + " \\\n  > " + normal_output),
                Job(mpileup_tumor_inputs, [tumor_output], command = "cat \\\n  " + " \\\n  ".join(mpileup_tumor_inputs) + " \\\n  > " + tumor_output)
            ], name = "rawmpileup_cat." + tumor_pair.name ))

        return jobs

    def paired_varscan2(self):
        """
        Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
        Varscan2 thresholds based on DREAM3 results generated by author see: https://github.com/dkoboldt/varscan/releases
        SSC INFO field remove to prevent collison with Samtools output during ensemble
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            for sequence in self.sequence_dictionary_variant():
                input_normal = os.path.join(varscan_directory, tumor_pair.normal.name + "." + sequence['name'] + ".mpileup")
                input_tumor = os.path.join(varscan_directory, tumor_pair.tumor.name + "." + sequence['name'] + ".mpileup")

                output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'])
                output_snp =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf")
                output_indel =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf")
                output_vcf = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf")
                output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf.gz")

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory], samples=[tumor_pair.normal, tumor_pair.tumor]),
                    varscan.somatic(input_normal, input_tumor, output, config.param('varscan2_somatic', 'other_options'), output_vcf_dep=output_vcf, output_snp_dep=output_snp, output_indel_dep=output_indel),
                    htslib.bgzip_tabix_vcf(output_snp, os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf.gz")),
                    htslib.bgzip_tabix_vcf(output_indel, os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf.gz")),
                    pipe_jobs([
                        bcftools.concat([os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf.gz"), os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf.gz")], None),
                        Job([None], [output_vcf], command="sed 's/TUMOR/"+ tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"+ tumor_pair.normal.name + "/g' | grep -v \"INFO=<ID=SSC\" | sed -E \"s/SSC=(.*);//g\" > " + output_vcf),
                    ]),
                    htslib.bgzip_tabix_vcf(output_vcf, output_vcf_gz),
                ], name = "varscan2_somatic." + tumor_pair.name + "." + sequence['name'] ) )

        return jobs

    def varscan2_fpfilter(self):
        """
        Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
        Varscan2 filtering thresholds based on DREAM3 results generated by author see: https://github.com/dkoboldt/varscan/releases
        Somatic and germline calls are filter using bam-readcount info from tumor and normal bam, respectively.
        Variants denoted as PASS and somatic (SS=1) or germline (SS=2) and loh (SS=3) are retained for further analysis.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            for sequence in self.sequence_dictionary_variant():
                input_vcf =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf")
                output_bed = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.bed")
                output_normal_readcount = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".normal.readcount")
                output_tumor_readcount =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".tumor.readcount")

                output_fpfilter_somatic =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.fpfilter.somatic.vcf")
                output_fpfilter_somatic_gz =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.fpfilter.somatic.vcf.gz")
                output_fpfilter_germline =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.fpfilter.germline_loh.vcf")
                output_fpfilter_germline_gz =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.fpfilter.germline_loh.vcf.gz")
                output_vcf_somatic_gz = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.somatic.vcf.gz")
                output_vcf_germline_loh_gz = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.germline_loh.vcf.gz")

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory,output_fpfilter_somatic,output_fpfilter_germline], samples=[tumor_pair.normal, tumor_pair.tumor]),
                    tools.vcf2bed(input_vcf, output_bed),
                    bamreadcount.readcount(inputNormal, output_bed, output_normal_readcount),
                    bamreadcount.readcount(inputTumor, output_bed, output_tumor_readcount),
                    varscan.fpfilter_somatic(input_vcf, output_tumor_readcount, output_fpfilter_somatic),
                    htslib.bgzip_tabix_vcf(output_fpfilter_somatic, output_fpfilter_somatic_gz),
                    pipe_jobs([
                        bcftools.view(output_fpfilter_somatic_gz, None, config.param('varscan2_readcount_fpfilter', 'somatic_filter_options')),
                        vcflib.vcfstreamsort(None, None),
                        htslib.bgzip_tabix_vcf(None, output_vcf_somatic_gz),
                    ]),
                    varscan.fpfilter_somatic(input_vcf, output_normal_readcount, output_fpfilter_germline),
                    htslib.bgzip_tabix_vcf(output_fpfilter_germline, output_fpfilter_germline_gz),
                    pipe_jobs([
                        bcftools.view(output_fpfilter_germline_gz, None, config.param('varscan2_readcount_fpfilter', 'germline_loh_filter_options')),
                        vcflib.vcfstreamsort(None, None),
                        htslib.bgzip_tabix_vcf(None, output_vcf_germline_loh_gz),
                    ])
                ], name = "varscan2_readcount_fpfilter." + tumor_pair.name + "." + sequence['name'] ) )

        return jobs

    def merge_varscan2(self):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            all_inputs = [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf.gz") for sequence in self.sequence_dictionary_variant()]

            #somatic_inputs = [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.somatic.vcf.gz") for sequence in self.sequence_dictionary_variant()]

            #germline_inputs = [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.germline_loh.vcf.gz") for sequence in self.sequence_dictionary_variant()]

            all_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz")
            all_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vt.vcf.gz")
            #somtic_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz")
            somtic_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            #germline_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vcf.gz")
            germline_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vt.vcf.gz")

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    bcftools.concat(all_inputs, None),
                    tools.fix_varscan_output(None, None),
                    Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                    Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                    Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),
                    htslib.bgzip_tabix_vcf(None, all_output),
                ]),
                pipe_jobs([
                    vt.decompose_and_normalize_mnps( all_output , None),
                    #vcflib.vcfstreamsort(None, None),
                    htslib.bgzip_tabix_vcf(None,  all_output_vt),
                ]),
                pipe_jobs([
                    bcftools.view(all_output_vt, None, config.param('varscan2_readcount_fpfilter', 'somatic_filter_options')),
                    htslib.bgzip_tabix_vcf(None, somtic_output_vt),
                ]),
                pipe_jobs([
                    bcftools.view(all_output_vt, None, config.param('varscan2_readcount_fpfilter', 'germline_loh_filter_options')),
                    htslib.bgzip_tabix_vcf(None, germline_output_vt),
                ]),
            ], name = "merge_varscan2." + tumor_pair.name ))

        return jobs


    def paired_mutect2(self):
        """
        GATK MuTect2 caller for SNVs and Indels.
        """

        jobs = []

        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")
            inputNormal = os.path.join("alignment",tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment",tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            mkdir_job = Job(command="mkdir -p " + mutect_directory, removable_files=[mutect_directory], samples=[tumor_pair.normal, tumor_pair.tumor])
            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    gatk.mutect2(inputNormal, inputTumor, os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz"))
                ], name="gatk_mutect2." + tumor_pair.name))

            else:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1)

                # Create one separate job for each of the first sequences
                for idx,sequences in enumerate(unique_sequences_per_job):
                    outPrefix =  tumor_pair.name + "." + str(idx) + ".mutect2"
                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        mkdir_job,
                        gatk.mutect2(inputNormal, inputTumor, os.path.join(mutect_directory, outPrefix + ".vcf.gz"), intervals=sequences)
                    ], name="gatk_mutect2." + tumor_pair.name + "." + str(idx)))

                # Create one last job to process the last remaining sequences and 'others' sequences
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    gatk.mutect2(inputNormal, inputTumor, os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"), exclude_intervals=unique_sequences_per_job_others)
                ], name="gatk_mutect2." + tumor_pair.name + ".others"))

        return jobs

    def merge_mutect2(self):
        """
        Merge SNVs and indels for mutect2
        Replace TUMOR and NORMAL sample names in vcf to the exact tumor/normal sample names
        Generate a somatic vcf containing only PASS variants
        """

        jobs = []

        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")
            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            output_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf.gz")
            output_vt_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vt.vcf.gz")
            output_somatic_vt = os.path.join(pair_directory, tumor_pair.name + ".mutect2.somatic.vt.vcf.gz")

            if nb_jobs == 1:
                input = os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz")
                jobs.append(concat_jobs([
                    Job([input], [output_gz], command="ln -s -f " + input + " " + output, samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(output_gz, None),
                        htslib.bgzip_tabix_vcf(None, output_vt_gz),
                    ]),
                    pipe_jobs([
                        Job([output_vt_gz], [None], command="zcat " + output + " | sed 's/TUMOR/"+ tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"+ tumor_pair.normal.name + "/g' | sed 's/Number=R/Number=./g' "),
                        bcftools.view(None, None, config.param('merge_filter_mutect2', 'filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_somatic),
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(output_somatic, None),
                        #vcflib.vcfstreamsort(None, None),
                        htslib.bgzip_tabix_vcf(None, output_somatic_vt),
                    ]),
                ], name="symlink_mutect_vcf." + tumor_pair.name))

            elif nb_jobs > 1:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1)

                # Create one separate job for each of the first sequences
                inputVCFs = []
                for idx,sequences in enumerate(unique_sequences_per_job):
                    inputVCFs.append(os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".mutect2.vcf.gz"))
                inputVCFs.append(os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"))

                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.concat(inputVCFs, None, config.param('merge_filter_mutect2', 'bcftools_options')),
                        Job([None], [None], command="sed 's/TUMOR/"+ tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"+ tumor_pair.normal.name + "/g' | sed 's/Number=R/Number=./g' " ),
                        htslib.bgzip_tabix_vcf(None, output_gz),
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(output_gz, None),
                        #vcflib.vcfstreamsort(None, None),
                        htslib.bgzip_tabix_vcf(None, output_vt_gz),
                    ]),
                    pipe_jobs([
                        bcftools.view(output_vt_gz, None, config.param('merge_filter_mutect2', 'filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_somatic_vt),
                     ]),
                 ], name = "merge_filter_mutect2." + tumor_pair.name))

        return jobs

    def samtools_paired(self):
        """
        Samtools caller for SNVs and Indels using verison 0.1.19.
        """

        jobs = []

        nb_jobs = config.param('samtools_paired', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            samtools_directory = os.path.join(pair_directory, "rawSamtools")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            paired_sample = [inputNormal,inputTumor]

            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + samtools_directory, removable_files=[samtools_directory], samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        samtools.mpileup(paired_sample, None, config.param('samtools_paired', 'mpileup_other_options'), ini_section="samtools_paired"),
                        samtools.bcftools_call_pair("-", os.path.join(samtools_directory,  tumor_pair.name + ".bcf"), config.param('samtools_paired', 'bcftools_view_options'), pair_calling=True),
                    ]),
                ], name="samtools_paired." + tumor_pair.name))

            else:
                for region in self.generate_approximate_windows(nb_jobs): #for idx,sequences in enumerate(unique_sequences_per_job):
                    jobs.append(concat_jobs([
                        Job(command="mkdir -p " + samtools_directory, removable_files=[samtools_directory], samples=[tumor_pair.normal, tumor_pair.tumor]),
                        pipe_jobs([
                            samtools.mpileup(paired_sample, None, config.param('samtools_paired', 'mpileup_other_options'), region,  ini_section="samtools_paired"),
                            samtools.bcftools_call_pair("-", os.path.join(samtools_directory,  tumor_pair.name + "." + region + ".bcf"), config.param('samtools_paired', 'bcftools_view_options'), pair_calling=True),
                        ]),
                ], name="samtools_paired." + tumor_pair.name + "." + region))

        return jobs

    def merge_filter_paired_samtools(self):
        """
        bcftools is used to merge the raw binary variants files created in the snpAndIndelBCF step.
        The output of bcftools is fed to varfilter, which does an additional filtering of the variants
        and transforms the output into the VCF (.vcf) format. One vcf file contain the SNP/INDEL calls
        for all samples in the experiment.
        Additional somatic filters are performed to reduce the number of FPs:
        1. vcflibs vcfsamplediff tags each variant with <tag>={germline,somatic,loh} to specify the type
        of variant given the genotype difference between the two samples.
        2. bcftools filter is used to retain only variants with CLR>=15 and have STATUS=somatic from
        vcfsamplediff
        3. bcftools filter is used to retain only variants that have STATUS=germline or STATUS=loh from
        vcfsamplediff
        """

        jobs = []
        nb_jobs = config.param('merge_filter_paired_samtools', 'approximate_nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            samtools_directory = os.path.join(pair_directory, "rawSamtools")
            output = os.path.join(samtools_directory,  tumor_pair.name + ".samtools.bcf")
            output_vcf = os.path.join(pair_directory,  tumor_pair.name + ".samtools.vcf.gz")
            output_vcf_vt = os.path.join(pair_directory,  tumor_pair.name + ".samtools.vt.vcf.gz")
            output_somatics = os.path.join(pair_directory,  tumor_pair.name + ".samtools.somatic.vt.vcf.gz")
            output_germline_loh = os.path.join(pair_directory,  tumor_pair.name + ".samtools.germline_loh.vt.vcf.gz")

            if nb_jobs == 1:
                inputs = os.path.join(samtools_directory,  tumor_pair.name + ".bcf")
                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    samtools.bcftools_cat_pair(inputs, output),
                    pipe_jobs([
                        samtools.bcftools_view_pair(output, None),
                        vcflib.vcfsamplediff(tumor_pair.normal.name, tumor_pair.tumor.name, None, None),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),
                        htslib.bgzip_tabix_vcf(None, output_vcf),
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(output_vcf, None),
                        #vcflib.vcfstreamsort(None, None),
                        htslib.bgzip_tabix_vcf(None, output_vcf_vt),
                    ]),
                    pipe_jobs([
                        bcftools.filter(output_vcf_vt, None, config.param('merge_filter_paired_samtools', 'somatic_filter_options')),
                        vcflib.vcffilter(None, None, config.param('merge_filter_paired_samtools', 'somatic_vcffilter_options')),
                        htslib.bgzip_tabix_vcf(None, output_somatics),
                    ]),
                    pipe_jobs([
                        bcftools.filter(output_vcf_vt, None, config.param('merge_filter_paired_samtools', 'germline_loh_filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_germline_loh),
                    ]),
                ], name = "merge_filter_paired_samtools." + tumor_pair.name))

            else:
                inputs = [os.path.join(samtools_directory,  tumor_pair.name + "." + region + ".bcf") for region in self.generate_approximate_windows(nb_jobs)] #for idx,sequences in enumerate(unique_sequences_per_job):
                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    samtools.bcftools_cat_pair(inputs, output),
                    pipe_jobs([
                        samtools.bcftools_view_pair(output, None),
                        vcflib.vcfsamplediff(tumor_pair.normal.name, tumor_pair.tumor.name, None, None),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),
                        htslib.bgzip_tabix_vcf(None, output_vcf),
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(output_vcf, None),
                        #vcflib.vcfstreamsort(None, None),
                        htslib.bgzip_tabix_vcf(None, output_vcf_vt),
                    ]),
                    pipe_jobs([
                        bcftools.filter(output_vcf_vt, None, config.param('merge_filter_paired_samtools', 'somatic_filter_options')),
                        #vcflib.vcffilter(None, None, config.param('merge_filter_paired_samtools', 'somatic_vcffilter_options')),
                        htslib.bgzip_tabix_vcf(None, output_somatics),
                    ]),
                     pipe_jobs([
                        bcftools.filter(output_vcf_vt, None, config.param('merge_filter_paired_samtools', 'germline_loh_filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_germline_loh),
                    ]),
                ], name = "merge_filter_paired_samtools." + tumor_pair.name))

        report_file = os.path.join("report", "DnaSeq.merge_filter_bcf.md")
        jobs.append(
            Job(
                [output_vcf],
                [report_file],
                command="""\
mkdir -p report && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_template_dir}/HumanVCFformatDescriptor.tsv \\
  report/ && \\
sed 's/\t/|/g' report/HumanVCFformatDescriptor.tsv | sed '2i-----|-----' >> {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="merge_filter_bcf_report")
        )

        return jobs

    def vardict_paired(self):
        """
        vardict caller for SNVs and Indels.
        Note: variants are filtered to remove instantance where REF == ALT and REF modified to 'N' when REF is AUPAC nomenclature
        """

        ##TO DO - the BED system needs to be revisted !!
        jobs = []

        nb_jobs = config.param('vardict_paired', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of vardict jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        use_bed = config.param('vardict_paired', 'use_bed', type='boolean', required=True)
        genome_dictionary = config.param('DEFAULT', 'genome_dictionary', type='filepath')
        if use_bed :
            bed = self.samples[0].readsets[0].beds[0]
            bed_intervals, interval_size = bed_file.parse_bed_file(bed)
            if not os.path.exists('vardict.tmp.0.bed'):
                bed_file_list = bed_file.split_by_size(bed_intervals, interval_size, nb_jobs, output="./vardict.tmp")
            else:
                bed_file_list = []
                for idx in range(nb_jobs):
                    bed_file_list.append(os.path.join("vardict.tmp." + str(idx) + ".bed"))

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            mkdir_job = Job(command="mkdir -p " + vardict_directory, removable_files=[vardict_directory], samples=[tumor_pair.normal, tumor_pair.tumor])

            if use_bed :
                idx = 0
                for bf in bed_file_list:
                    output=os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz")
                    jobs.append(concat_jobs([
                        mkdir_job,
                        pipe_jobs([
                            vardict.paired_java(inputNormal, inputTumor, tumor_pair.name, None, bf),
                            vardict.testsomatic(None,None),
                            vardict.var2vcf(None, tumor_pair.normal.name, tumor_pair.tumor.name, None ),
                            htslib.bgzip_tabix_vcf(None, output),
                        ]),
                    ], name="vardict_paired." + tumor_pair.name + "." + str(idx) ))
                    idx += 1
            else:
                beds = []
                for idx in range(nb_jobs):
                    beds.append(os.path.join(vardict_directory, "chr." + str(idx) + ".bed"))
                if nb_jobs == 1:
                    bedJob = vardict.dict2beds(genome_dictionary, beds)
                    output =os.path.join(vardict_directory, tumor_pair.name + ".0.vardict.vcf.gz")
                    jobs.append(concat_jobs([
                        mkdir_job,
                        bedJob,
                        pipe_jobs([
                            vardict.paired_java(inputNormal, inputTumor, tumor_pair.name, None, beds.pop()),
                            vardict.testsomatic(None,None),
                            vardict.var2vcf(None, tumor_pair.normal.name, tumor_pair.tumor.name, None ),
                            htslib.bgzip_tabix_vcf(None, output),
                        ]),
                     ], name="vardict_paired." + tumor_pair.name + ".0"))
                else:
                    bedJob = vardict.dict2beds(genome_dictionary, beds)
                    jobs.append(concat_jobs([
                        mkdir_job,
                        bedJob
                    ], name="vardict.genome.beds." + tumor_pair.name))

                    for idx in range(nb_jobs):
                        output = os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz")
                        jobs.append(concat_jobs([
                            mkdir_job,
                            pipe_jobs([
                                vardict.paired_java(inputNormal, inputTumor, tumor_pair.name, None, beds[idx]),
                                vardict.testsomatic(None, None),
                                vardict.var2vcf(None, tumor_pair.normal.name, tumor_pair.tumor.name, None),
                                htslib.bgzip_tabix_vcf(None, output),
                            ]),
                        ], name="vardict_paired." + tumor_pair.name + "." + str(idx) ))
        return jobs

    def merge_filter_paired_vardict(self):
        """
        The fully merged vcf is filtered using following steps:
        1. Retain only variants designated as somatic by VarDict: either StrongSomatic or LikelySomatic
        2. Somatics identified in step 1 must have PASS filter
        """

        jobs = []
        nb_jobs = config.param('vardict_paired', 'nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            output_tmp = os.path.join(pair_directory,  tumor_pair.name + ".vardict.tmp.vcf.gz")
            output = os.path.join(pair_directory,  tumor_pair.name + ".vardict.vcf.gz")
            output_vt = os.path.join(pair_directory,  tumor_pair.name + ".vardict.vt.vcf.gz")
            output_somatic = os.path.join(pair_directory,  tumor_pair.name + ".vardict.somatic.vt.vcf.gz")
            output_germline_loh = os.path.join(pair_directory,  tumor_pair.name + ".vardict.germline_loh.vt.vcf.gz")

            if nb_jobs == 1:
                inputs = os.path.join(vardict_directory,  tumor_pair.name + ".0.vardict.vcf.gz")
                jobs.append(concat_jobs([
                    Job([input], [output_tmp], command="ln -s -f " + input + " " + output_tmp, samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        Job([output_tmp], [None], command="zcat {output} | awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),
                        htslib.bgzip_tabix_vcf(None, output)
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(output, None),
                        vcflib.vcfstreamsort(None, None),
                        htslib.bgzip_tabix_vcf(None, output_vt),
                    ]),
                    pipe_jobs([
                       bcftools.view(output_vt, None, config.param('merge_filter_paired_vardict', 'somatic_filter_options')),
                       htslib.bgzip_tabix_vcf(None, output_somatics),
                    ]),
                    pipe_jobs([
                       bcftools.view(output_vt, None, config.param('merge_filter_paired_vardict', 'germline_loh_filter_options')),
                       htslib.bgzip_tabix_vcf(None, output_germline_loh),
                    ]),
                ], removable_files=[output], name="symlink_vardict_vcf." + tumor_pair.name ))
            else:
                inputVCFs = []
                for idx in range(nb_jobs):
                    inputVCFs.append(os.path.join(vardict_directory,  tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz"))

                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.concat(inputVCFs, None),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),
                        htslib.bgzip_tabix_vcf(None, output),
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(output, None),
                        vcflib.vcfstreamsort(None, None),
                        htslib.bgzip_tabix_vcf(None, output_vt),
                    ]),
                    pipe_jobs([
                       bcftools.view(output_vt, None, config.param('merge_filter_paired_vardict', 'somatic_filter_options')),
                       htslib.bgzip_tabix_vcf(None, output_somatic),
                    ]),
                    pipe_jobs([
                       bcftools.view(output_vt, None, config.param('merge_filter_paired_vardict', 'germline_loh_filter_options')),
                       htslib.bgzip_tabix_vcf(None, output_germline_loh),
                    ]),
                ], name = "merge_filter_paired_vardict." + tumor_pair.name))

        return jobs


    def ensemble_somatic(self):
        """
        Apply Bcbio.variations ensemble approach for mutect2, Vardict, Samtools and VarScan2 calls
        Filter ensemble calls to retain only calls overlapping 2 or more callers
        """

        jobs = []
        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_directory = os.path.join("pairedVariants", tumor_pair.name)

            input_mutect2 = os.path.join(input_directory, tumor_pair.name + ".mutect2.vt.vcf.gz")
            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.somatic.vt.vcf.gz")
            input_samtools = os.path.join(input_directory, tumor_pair.name + ".samtools.somatic.vt.vcf.gz")
            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            inputSomaticVCFs = [input_mutect2, input_vardict, input_samtools, input_varscan2]

            output_ensemble = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.vt.vcf.gz")
            output_flt = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.flt.vcf.gz")

            mkdir_job = Job(command="mkdir -p " + paired_ensemble_directory, removable_files=[os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic-work"), output_ensemble], samples=[tumor_pair.normal, tumor_pair.tumor])

            jobs.append(concat_jobs([
                # Create output directory since it is not done by default by GATK tools
                mkdir_job,
                bcbio_variation_recall.ensemble(inputSomaticVCFs, output_ensemble, config.param('bcbio_ensemble_somatic', 'options')),
            ], name="bcbio_ensemble_somatic." + tumor_pair.name))

        return jobs


    def ensemble_germline_loh(self):
        """
        Apply Bcbio.variations ensemble approach for Vardict, Samtools and VarScan2 calls
        Filter ensemble calls to retain only calls overlapping 2 or more callers
        """

        jobs = []
        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_directory = os.path.join("pairedVariants", tumor_pair.name)

            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.germline_loh.vt.vcf.gz")
            input_samtools = os.path.join(input_directory, tumor_pair.name + ".samtools.germline_loh.vt.vcf.gz")
            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.germline_loh.vt.vcf.gz")
            inputGermlineVCFs = [input_vardict, input_samtools, input_varscan2]

            output_ensemble = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline_loh.vt.vcf.gz")

            mkdir_job = Job(command="mkdir -p " + paired_ensemble_directory, removable_files=[os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline_loh-work"), output_ensemble], samples=[tumor_pair.normal, tumor_pair.tumor])

            jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    bcbio_variation_recall.ensemble(inputGermlineVCFs, output_ensemble, config.param('bcbio_ensemble_germline_loh', 'options')),
                ], name="bcbio_ensemble_germline_loh." + tumor_pair.name))

        return jobs

    def gatk_variant_annotator_somatic(self):
        """
        Add vcf annotations to ensemble vcf: Standard and Somatic annotations
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            inputNormal = os.path.join("alignment",tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment",tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            input_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.vcf.gz")

            for sequence in self.sequence_dictionary_variant():
                output_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot." + sequence['name'] + ".vcf.gz")

                mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output_somatic_variants], samples=[tumor_pair.normal, tumor_pair.tumor])

                jobs.append(concat_jobs([
                    mkdir_job,
                    gatk.variant_annotator( inputNormal, inputTumor, input_somatic_variants, output_somatic_variants, intervals=[sequence['name']] ),
                ], name = "gatk_variant_annotator.somatic." + sequence['name'] + "." + tumor_pair.name))

        return jobs

    def gatk_variant_annotator_germline(self):
        """
        Add vcf annotations to ensemble vcf: most importantly the AD field
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            inputNormal = os.path.join("alignment",tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment",tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            input_germline_loh_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline_loh.vt.vcf.gz")

            for sequence in self.sequence_dictionary_variant():
                output_germline_loh_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline_loh.vt.annot." + sequence['name'] + ".vcf.gz")

                mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output_germline_loh_variants], samples=[tumor_pair.normal, tumor_pair.tumor])

                jobs.append(concat_jobs([
                    mkdir_job,
                    gatk.variant_annotator(inputNormal, inputTumor, input_germline_loh_variants, output_germline_loh_variants, intervals=[sequence['name']]),
                ], name = "gatk_variant_annotator.germline." + sequence['name'] + "." + tumor_pair.name))

        return jobs

    def merge_gatk_variant_annotator_somatic(self):
        """
        Merge annotated somatic vcfs
        """
        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():

            somatic_inputs = [os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot." + sequence['name'] + ".vcf.gz") for sequence in self.sequence_dictionary_variant()]
            somatic_output = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    bcftools.concat(somatic_inputs, None),
                    htslib.bgzip_tabix_vcf(None, somatic_output),
                ]),
            ], name = "merge_gatk_variant_annotator.somatic." + tumor_pair.name))

        return jobs

    def merge_gatk_variant_annotator_germline(self):
        """
        Merge annotated germline and LOH vcfs
        """
        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():

            germline_inputs = [os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline_loh.vt.annot." + sequence['name'] + ".vcf.gz") for sequence in self.sequence_dictionary_variant()]
            germline_output = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline_loh.vt.annot.vcf.gz")

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    bcftools.concat(germline_inputs, None),
                    htslib.bgzip_tabix_vcf(None, germline_output),
                ]),
            ], name = "merge_gatk_variant_annotator.germline." + tumor_pair.name))

        return jobs

    def compute_cancer_effects_somatic(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        if not os.path.exists(ensemble_directory):
            os.makedirs(ensemble_directory)

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            if not os.path.exists(paired_directory):
                os.makedirs(paired_directory)

            input_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")
            output_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.vcf")
            output_somatic_gz = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.vcf.gz")

            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write( tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n" )

            mkdir_job = Job(command="mkdir -p " + paired_directory, removable_files=[output_somatic], samples=[tumor_pair.normal, tumor_pair.tumor])

            jobs.append(concat_jobs([
                mkdir_job,
                snpeff.compute_effects(input_somatic, output_somatic, cancer_sample_file=cancer_pair_filename, options=config.param('compute_cancer_effects_somatic', 'options')),
                htslib.bgzip_tabix_vcf(output_somatic, output_somatic_gz),
            ], name = "compute_cancer_effects_somatic." + tumor_pair.name))

        return jobs

    def compute_cancer_effects_germline(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)

            input_germline = os.path.join(paired_directory, tumor_pair.name + ".ensemble.germline_loh.vt.annot.vcf.gz")
            output_germline = os.path.join(paired_directory, tumor_pair.name + ".ensemble.germline_loh.vt.annot.snpeff.vcf")
            output_germline_gz = os.path.join(paired_directory, tumor_pair.name + ".ensemble.germline_loh.vt.annot.snpeff.vcf.gz")

            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write( tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n" )

            mkdir_job = Job(command="mkdir -p " + paired_directory, removable_files=[output_germline], samples=[tumor_pair.normal, tumor_pair.tumor])

            jobs.append(concat_jobs([
                mkdir_job,
                snpeff.compute_effects(input_germline, output_germline, options=config.param('compute_cancer_effects_germline', 'options')),
                htslib.bgzip_tabix_vcf(output_germline, output_germline_gz),
            ], name = "compute_cancer_effects_germline." + tumor_pair.name))

        return jobs

    def combine_tumor_pairs_somatic(self):
        """
        Combine numerous ensemble vcfs into one vcf for gemini annotations
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input_merged_vcfs = [os.path.join(ensemble_directory, tumor_pair.name , tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz") for tumor_pair in self.tumor_pairs.itervalues()]
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output], samples=self.samples)

        if len(input_merged_vcfs) == 1:
            jobs.append(concat_jobs([
                mkdir_job,
                Job([input_merged_vcfs[0]], [output], command="ln -s -f " + os.path.abspath(input_merged_vcfs[0]) + " " + output)
            ],name="gatk_combine_variants.somatic.allPairs"))

        else:
            jobs.append(concat_jobs([
                mkdir_job,
                gatk.combine_variants(input_merged_vcfs, output)
            ], name="gatk_combine_variants.somatic.allPairs"))

        return jobs

    def combine_tumor_pairs_germline(self):
        """
        Combine numerous ensemble vcfs into one vcf for gemini annotations
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input_merged_vcfs = [os.path.join(ensemble_directory, tumor_pair.name , tumor_pair.name + ".ensemble.germline_loh.vt.annot.vcf.gz") for tumor_pair in self.tumor_pairs.itervalues()]
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.vcf.gz")

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output], samples=self.samples)

        if len(input_merged_vcfs) == 1:
            jobs.append(concat_jobs([
                mkdir_job,
                Job([input_merged_vcfs[0]], [output], command="ln -s -f " + os.path.abspath(input_merged_vcfs[0]) + " " + output)
            ],name="gatk_combine_variants.germline_loh.allPairs"))

        else:

            jobs.append(concat_jobs([
                mkdir_job,
                gatk.combine_variants(input_merged_vcfs, output)
            ], name="gatk_combine_variants.germline_loh.allPairs"))

        return jobs

    def decompose_and_normalize_mnps_somatic(self):
        """
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt)
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output], samples=self.samples)

        jobs.append(concat_jobs([
            mkdir_job,
            vt.decompose_and_normalize_mnps(input, output)
        ],name = "decompose_and_normalize_mnps.somatic.allPairs"))

        return jobs

    def decompose_and_normalize_mnps_germline(self):
        """
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt)
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.vcf.gz")

        job = vt.decompose_and_normalize_mnps(input, output)
        job.name = "decompose_and_normalize_mnps.germline.allPairs"

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output], samples=self.samples)

        jobs.append(concat_jobs([
            mkdir_job,
            vt.decompose_and_normalize_mnps(input, output)
        ],name = "decompose_and_normalize_mnps.somatic.allPairs"))

        return jobs

    def all_pairs_compute_effects_somatic(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        Applied to all tumor pairs.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input =  os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.snpeff.vcf")
        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.snpeff.vcf.gz")

        cancer_pair_filename = os.path.join('cancer_snpeff.tsv')
        cancer_pair = open(cancer_pair_filename, 'w')

        for tumor_pair in self.tumor_pairs.itervalues():
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output], samples=self.samples)

        jobs.append(concat_jobs([
            mkdir_job,
            snpeff.compute_effects(input, output, cancer_sample_file=cancer_pair_filename,options=config.param('compute_cancer_effects_somatic', 'options')),
            htslib.bgzip_tabix_vcf(output, output_gz),
        ],name = "compute_effects.somatic.allPairs"))

        return jobs

    def all_pairs_compute_effects_germline(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        Applied to all tumor pairs.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input =  os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.snpeff.vcf")
        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.snpeff.vcf.gz")

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output], samples=self.samples)

        jobs.append(concat_jobs([
            mkdir_job,
            snpeff.compute_effects(input, output, options=config.param('compute_cancer_effects_germline', 'options')),
            htslib.bgzip_tabix_vcf(output, output_gz),
        ], name = "compute_effects.germline.allPair"))

        return jobs

    def gemini_annotations_somatic(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_prefix = os.path.join(ensemble_directory, "allPairs")
        gemini_module=config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2],gemini_module[-1]])

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + ensemble_directory, samples=self.samples),
            gemini.gemini_annotations( gemini_prefix + ".ensemble.somatic.vt.annot.snpeff.vcf.gz", gemini_prefix + ".somatic.gemini." + gemini_version + ".db", temp_dir)
        ], name="gemini_annotations.somatic.allPairs"))

        return jobs

    def gemini_annotations_germline(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_prefix = os.path.join(ensemble_directory, "allPairs")
        gemini_module=config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2],gemini_module[-1]])

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + ensemble_directory, samples=self.samples),
            gemini.gemini_annotations( gemini_prefix + ".ensemble.germline_loh.vt.annot.snpeff.vcf.gz", gemini_prefix + ".germline_loh.gemini." + gemini_version + ".db", temp_dir)
        ], name="gemini_annotations.germline.allPairs"))

        return jobs

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.bwa_mem_picard_sort_sam,
            self.sambamba_merge_sam_files,
            self.gatk_indel_realigner,
            self.sambamba_merge_realigned,
            #self.fix_mate_by_coordinate,
            self.sambamba_mark_duplicates,
            self.recalibration,
            self.conpair_concordance_contamination,
            self.rawmpileup_panel,
            self.paired_varscan2_panel,
            self.merge_varscan2_panel,
            self.preprocess_vcf_panel,
            self.snp_effect_panel,
            self.gemini_annotations_panel,
            self.metrics,
            self.picard_calculate_hs_metrics,
            self.gatk_callable_loci,
            self.extract_common_snp_freq,
            self.baf_plot,
            self.rawmpileup,
            #self.rawmpileup_cat,
            self.paired_varscan2,
            #self.varscan2_fpfilter,
            self.merge_varscan2,
            self.paired_mutect2,
            self.merge_mutect2,
            self.samtools_paired,
            self.merge_filter_paired_samtools,
            self.vardict_paired,
            self.merge_filter_paired_vardict,
            self.ensemble_somatic,
            self.gatk_variant_annotator_somatic,
            self.merge_gatk_variant_annotator_somatic,
            self.compute_cancer_effects_somatic,
            self.combine_tumor_pairs_somatic,
            #self.decompose_and_normalize_mnps_somatic,
            self.all_pairs_compute_effects_somatic,
            self.gemini_annotations_somatic,
            self.ensemble_germline_loh,
            self.gatk_variant_annotator_germline,
            self.merge_gatk_variant_annotator_germline,
            self.compute_cancer_effects_germline,
            self.combine_tumor_pairs_germline,
            #self.decompose_and_normalize_mnps_germline,
            self.all_pairs_compute_effects_germline,
            self.gemini_annotations_germline
        ]

if __name__ == '__main__':
    TumorPair()
