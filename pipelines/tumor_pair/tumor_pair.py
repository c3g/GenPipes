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

from bfx import bcftools
from bfx import gq_seq_utils
from bfx import gatk
from bfx import picard
from bfx import samtools
from bfx import scalpel
from bfx import tools
from pipelines.dnaseq import dnaseq

log = logging.getLogger(__name__)

class TumorPair(dnaseq.DnaSeq):
    """
    Tumor Pair Pipeline
    =================

    The Tumor Pair pipeline is based on the DNA-Seq pipeline and follow the first initial steps.
    The difference starts at the variant calling step. At this point, a paired caller is used to call SNVs and Indels
    from the pairs given as input.
    """

    def __init__(self):
        # Add pipeline specific arguments
        self.argparser.add_argument("-p", "--pairs", help="pairs file", type=file)
        super(TumorPair, self).__init__()

    @property
    def tumor_pairs(self):
        if not hasattr(self, "_tumor_pairs"):
            self._tumor_pairs = parse_tumor_pair_file(self.args.pairs.name, self.samples)
        return self._tumor_pairs

    def paired_SNVs(self):
        """
        GATK MuTect caller for SNVs.
        """

        jobs = []

        nb_jobs = config.param('gatk_mutect', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect")
            inputNormal = os.path.join("alignment",tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment",tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            mkdir_job = Job(command="mkdir -p " + mutect_directory)
            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    gatk.mutect(inputNormal, inputTumor, os.path.join(mutect_directory, tumor_pair.name + ".mutect.stats"), os.path.join(mutect_directory, tumor_pair.name + ".mutect.vcf"))
                ], name="gatk_mutect." + tumor_pair.name))

            else:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                # Create one separate job for each of the first sequences
                for idx,sequences in enumerate(unique_sequences_per_job):
                    outPrefix =  tumor_pair.name + "." + str(idx) + ".mutect"
                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        mkdir_job,
                        gatk.mutect(inputNormal, inputTumor, os.path.join(mutect_directory, outPrefix + ".stats"), os.path.join(mutect_directory, outPrefix + ".vcf"), intervals=sequences)
                    ], name="gatk_mutect." + tumor_pair.name + "." + str(idx)))

                # Create one last job to process the last remaining sequences and 'others' sequences
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    gatk.mutect(inputNormal, inputTumor, os.path.join(mutect_directory, tumor_pair.name + ".others.mutect.stats"), os.path.join(mutect_directory, tumor_pair.name + ".others.mutect.vcf"), exclude_intervals=unique_sequences_per_job_others)
                ], name="gatk_mutect." + tumor_pair.name + ".others"))

        return jobs

    def merge_SNVs(self):
        """
        Merge SNVs
        """

        jobs = []

        nb_jobs = config.param('gatk_mutect', 'nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect")
            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            output = os.path.join(pair_directory, tumor_pair.name + ".mutect.vcf")
            if nb_jobs == 1:
                input = os.path.join(mutect_directory, tumor_pair.name + ".mutect.vcf")
                job = Job([input], [output], command="ln -s -f " + input + " " + output, removable_files=[output], name="symlink_mutect_vcf." + tumor_pair.name)

            elif nb_jobs > 1:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                # Create one separate job for each of the first sequences
                inputVCFs = []
                for idx,sequences in enumerate(unique_sequences_per_job):
                    inputVCFs.append(os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".mutect.vcf"))
                inputVCFs.append(os.path.join(mutect_directory, tumor_pair.name + ".others.mutect.vcf"))

                job = gatk.cat_variants(inputVCFs, output)
                job.name = "gatk_cat_mutect." + tumor_pair.name

            jobs.append(job)

        return jobs

    def paired_indels(self):
        """
        Scalpel caller for insertions and deletions.
        """

        jobs = []

        nb_jobs = config.param('scalpel', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of scalpel jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            scalpel_directory = os.path.join(pair_directory, "rawScalpel")
            inputNormal = os.path.join("alignment",tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment",tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            mkdir_job = Job(command="mkdir -p " + scalpel_directory)

            genome_dictionary = config.param('DEFAULT', 'genome_dictionary', type='filepath')
            beds = []
            for idx in range(nb_jobs):
                beds.append(os.path.join(scalpel_directory, 'chrs.' + str(idx) + '.bed'))
            if nb_jobs == 1:
                bedJob = tools.dict2beds(genome_dictionary, beds)
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    bedJob,
                    scalpel.scalpel_somatic(inputNormal, inputTumor, os.path.join(scalpel_directory, tumor_pair.name + ".scalpel"), beds.pop())
                ], name="scalpel." + tumor_pair.name))

            else:
                bedJob = tools.dict2beds(genome_dictionary, beds)
                jobs.append(concat_jobs([mkdir_job,bedJob], name="scalpel.genome.beds." + tumor_pair.name))

                for idx in range(nb_jobs):
                    scalpelJob = scalpel.scalpel_somatic(inputNormal, inputTumor, os.path.join(scalpel_directory, tumor_pair.name + '.' + str(idx) + '.scalpel'), beds[idx])
                    scalpelJob.name = "scalpel." + tumor_pair.name + "." + str(idx)
                    jobs.append(scalpelJob)
        return jobs

    def merge_indels(self):
        """
        Merge Insertion and Deletions
        """

        jobs = []

        nb_jobs = config.param('scalpel', 'nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            scalpel_directory = os.path.join(pair_directory, "rawScalpel")
            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            output = os.path.join(pair_directory, tumor_pair.name + ".scalpel.vcf")
            if nb_jobs == 1:
                outputWrongName = os.path.join(scalpel_directory, tumor_pair.name + ".rename.scalpel.vcf")
                outputSomatic = os.path.join(scalpel_directory, tumor_pair.name + ".scalpel.somatic.vcf")
                outputCommon = os.path.join(scalpel_directory, tumor_pair.name + ".scalpel.common.vcf")
                input_somatic = os.path.join(scalpel_directory, tumor_pair.name + ".scalpel", 'main', 'somatic.5x.indel.vcf')
                input_common = os.path.join(scalpel_directory, tumor_pair.name + ".scalpel", 'main', 'common.5x.indel.vcf')

                job = concat_jobs([
                    bcftools.add_reject(input_common, outputCommon),
                    bcftools.add_chi2(input_somatic, outputSomatic),
                    picard.sort_vcfs([outputCommon, outputSomatic], outputWrongName),
                    Job([outputWrongName], [output], command="sed 's/sample_name/" + tumor_pair.tumor.name + "/g'" + outputWrongName + ' > ' + output)]
                    , removable_files=[outputCommon,outputSomatic,outputWrongName,output]
                    , name="scalpel_merged_vcf." + tumor_pair.name)
                jobs.append(job)

            elif nb_jobs > 1:
                vcfsToMerge =[]
                for idx in range(nb_jobs):
                    sample_directory = os.path.join(scalpel_directory, tumor_pair.name + "." + str(idx) + ".scalpel")
                    outputWrongName = os.path.join(sample_directory, tumor_pair.name + ".rename.scalpel.vcf")
                    outputSomatic = os.path.join(sample_directory, tumor_pair.name + ".scalpel.somatic.vcf")
                    outputCommon = os.path.join(sample_directory, tumor_pair.name + ".scalpel.common.vcf")
                    outputIdx = os.path.join(sample_directory, tumor_pair.name + ".scalpel.vcf")
                    input_somatic = os.path.join(sample_directory, 'main', 'somatic.5x.indel.vcf')
                    input_common = os.path.join(sample_directory, 'main', 'common.5x.indel.vcf')
                    
                    vcfsToMerge.append(outputIdx)
                    job = concat_jobs([
                        bcftools.add_reject(input_common, outputCommon),
                        bcftools.add_chi2Filter(input_somatic, outputSomatic),
                        picard.sort_vcfs([outputCommon, outputSomatic], outputWrongName),
                        Job([outputWrongName], [output], command="sed 's/sample_name/"+ tumor_pair.tumor.name + "/g'" + outputWrongName + ' > ' + outputIdx, removable_files=[outputCommon,outputSomatic,outputWrongName])]
                        , name="scalpel_merged_vcf." + tumor_pair.name + "." + str(idx))
                    jobs.append(job)
                job = gatk.cat_variants(vcfsToMerge, output)
                job.name="scalpel_merge_all_vcfs." + tumor_pair.name
                jobs.append(job)

        return jobs

    def merge_svns_and_indels(self):
        """
        Merge snvs and indels
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            input_mutect = os.path.join(pair_directory, tumor_pair.name + ".mutect.vcf")
            input_scalpel = os.path.join(pair_directory, tumor_pair.name + ".scalpel.vcf")
            output = os.path.join(pair_directory, tumor_pair.name + ".merged.vcf")

            job = concat_jobs([
                gatk.combine_variants([input_mutect, input_scalpel], output)],
                name="merge_snvs_and_indels." + tumor_pair.name)
            jobs.append(job)

        return jobs

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
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
            self.paired_SNVs,
            self.merge_SNVs,
            self.paired_indels,
            self.merge_indels,
            self.merge_svns_and_indels
        ]

if __name__ == '__main__': 
    TumorPair()
