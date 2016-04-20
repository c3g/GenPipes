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
from bfx import bed_file
from bfx import vardict
from bfx import bcbio_variation
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

            mkdir_job = Job(command="mkdir -p " + mutect_directory)
            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    gatk.mutect2(inputNormal, inputTumor, os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf"))
                ], name="gatk_mutect2." + tumor_pair.name))

            else:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                # Create one separate job for each of the first sequences
                for idx,sequences in enumerate(unique_sequences_per_job):
                    outPrefix =  tumor_pair.name + "." + str(idx) + ".mutect2"
                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        mkdir_job,
                        gatk.mutect2(inputNormal, inputTumor, os.path.join(mutect_directory, outPrefix + ".vcf"), intervals=sequences)
                    ], name="gatk_mutect2." + tumor_pair.name + "." + str(idx)))

                # Create one last job to process the last remaining sequences and 'others' sequences
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    gatk.mutect2(inputNormal, inputTumor, os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf"), exclude_intervals=unique_sequences_per_job_others)
                ], name="gatk_mutect2." + tumor_pair.name + ".others"))

        return jobs

    def merge_mutect2(self):
        """
        Merge SNVs and indels for mutect2
        """

        jobs = []

        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")
            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            output = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf")
            if nb_jobs == 1:
                input = os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf")
                job = Job([input], [output], command="ln -s -f " + input + " " + output, removable_files=[output], name="symlink_mutect_vcf." + tumor_pair.name)

            elif nb_jobs > 1:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                # Create one separate job for each of the first sequences
                inputVCFs = []
                for idx,sequences in enumerate(unique_sequences_per_job):
                    inputVCFs.append(os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".mutect2.vcf"))
                inputVCFs.append(os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf"))

                job = gatk.cat_variants(inputVCFs, output)
                job.name = "gatk_cat_mutect." + tumor_pair.name

            jobs.append(job)

        return jobs

    def samtools_paired(self):
        """
        Samtools caller for SNVs and Indels.
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
                    Job(command="mkdir -p " + samtools_directory),
                    pipe_jobs([
                        samtools.mpileup(paired_sample, None, config.param('samtools_paired', 'mpileup_other_options')),
                        samtools.bcftools_view("-", os.path.join(samtools_directory,  tumor_pair.name + ".bcf"), config.param('samtools_paired', 'bcftools_view_options'), pair_calling=True),
                    ])], name="samtools_paired." + tumor_pair.name))

            else:
                for region in self.generate_approximate_windows(nb_jobs): #for idx,sequences in enumerate(unique_sequences_per_job):
                    jobs.append(concat_jobs([
                        Job(command="mkdir -p " + samtools_directory),
                        pipe_jobs([
                            samtools.mpileup(paired_sample, None, config.param('samtools_paired', 'mpileup_other_options'), region),
                            samtools.bcftools_view("-", os.path.join(samtools_directory,  tumor_pair.name + "." + region + ".bcf"), config.param('samtools_paired', 'bcftools_view_options'), pair_calling=True),
                        ])], name="samtools_paired." + tumor_pair.name + "." + region))

        return jobs

    def merge_filter_paired_samtools(self):
        """
        bcftools is used to merge the raw binary variants files created in the snpAndIndelBCF step.
        The output of bcftools is fed to varfilter, which does an additional filtering of the variants
        and transforms the output into the VCF (.vcf) format. One vcf file contain the SNP/INDEL calls
        for all samples in the experiment.
        """

        jobs = []
        nb_jobs = config.param('merge_filter_paired_samtools', 'approximate_nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            samtools_directory = os.path.join(pair_directory, "rawSamtools")
            output = os.path.join(samtools_directory,  tumor_pair.name + ".samtools.bcf")
            output_vcf = os.path.join(pair_directory,  tumor_pair.name + ".samtools.vcf")

            if nb_jobs == 1:
                inputs = os.path.join(samtools_directory,  tumor_pair.name + ".bcf") 
                jobs.append(concat_jobs([
                    samtools.bcftools_cat(inputs, output),
                    samtools.bcftools_view(output, output_vcf)
                ], name = "merge_filter_paired_samtools." + tumor_pair.name))

            else:
                inputs = [os.path.join(samtools_directory,  tumor_pair.name + "." + region + ".bcf") for region in self.generate_approximate_windows(nb_jobs)] #for idx,sequences in enumerate(unique_sequences_per_job):
                jobs.append(concat_jobs([
                    samtools.bcftools_cat(inputs, output),
                    samtools.bcftools_view(output, output_vcf)
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
        """

        jobs = []

        nb_jobs = config.param('vardict_paired', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of vardict jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        use_bed = config.param('vardict_paired', 'use_bed', type='boolean', required=True)
        genome_dictionary = config.param('DEFAULT', 'genome_dictionary', type='filepath')
        if use_bed :
            bed = self.samples[0].readsets[0].beds[0]
            bed_intervals, interval_size = bed_file.parse_bed_file(bed)
            bed_file_list = bed_file.split_by_size(bed_intervals, interval_size, nb_jobs, output="pairedVariants/vardict.tmp")

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            mkdir_job = Job(command="mkdir -p " + vardict_directory)

            if use_bed :
                idx = 0
                for bf in bed_file_list:
                    jobs.append(concat_jobs([
                        mkdir_job,
                        pipe_jobs([
                            vardict.paired_java(inputNormal, inputTumor, tumor_pair.name, None, bf),
                            vardict.testsomatic(None,None),
                            vardict.var2vcf(os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf"), tumor_pair.normal.name, tumor_pair.tumor.name, None ),
                    ])], name="vardict_paired." + tumor_pair.name + "." + str(idx) ))
                    idx += 1
            else:
                beds = []
                for idx in range(nb_jobs):
                    beds.append(os.path.join(vardict_directory, "chr." + str(idx) + ".bed"))
                if nb_jobs == 1:
                    bedJob = vardict.dict2beds(genome_dictionary, beds)
                    jobs.append(concat_jobs([
                        mkdir_job,
                        bedJob,
                        pipe_jobs([
                            vardict.paired_java(inputNormal, inputTumor, tumor_pair.name, None, beds.pop()),
                            vardict.testsomatic(None,None),
                            vardict.var2vcf(os.path.join(vardict_directory, tumor_pair.name + ".0.vardict.vcf"), tumor_pair.normal.name, tumor_pair.tumor.name, None ),
                        ])], name="vardict_paired." + tumor_pair.name + ".0"))
                else:
                    bedJob = vardict.dict2beds(genome_dictionary, beds)
                    jobs.append(concat_jobs([mkdir_job,bedJob], name="vardict.genome.beds." + tumor_pair.name))

                    for idx in range(nb_jobs):
                        jobs.append(concat_jobs([
                        mkdir_job,
                        pipe_jobs([
                            vardict.paired_java(inputNormal, inputTumor, tumor_pair.name, None, beds[idx]),
                            vardict.testsomatic(None, None),
                            vardict.var2vcf(os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf"), tumor_pair.normal.name, tumor_pair.tumor.name, None),
                        ])], name="vardict_paired." + tumor_pair.name + "." + str(idx) ))                        
        return jobs

    def merge_paired_vardict(self):
        """
        
        """

        jobs = []
        nb_jobs = config.param('vardict_paired', 'nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            output = os.path.join(pair_directory,  tumor_pair.name + ".vardict.vcf")

            if nb_jobs == 1:
                inputs = os.path.join(vardict_directory,  tumor_pair.name + ".0.vardict.vcf")
                job = Job([input], [output], command="ln -s -f " + input + " " + output, removable_files=[output], name="symlink_vardict_vcf." + tumor_pair.name)
                
            else:
                inputVCFs = []
                for idx in range(nb_jobs):
                    inputVCFs.append(os.path.join(vardict_directory,  tumor_pair.name + "." + str(idx) + ".vardict.vcf"))
                
                job = gatk.cat_variants(inputVCFs, output)
                job.name = "gatk_cat_vardict." + tumor_pair.name

        jobs.append(job)

        return jobs


    def ensemble_tumor_pairs(self):
        """
        """

        jobs = []
        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            input_directory = os.path.join("pairedVariants", tumor_pair.name)
            input_mutect2 = os.path.join(input_directory, tumor_pair.name + ".mutect2.vcf")
            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.vcf")
            input_samtools = os.path.join(input_directory, tumor_pair.name + ".samtools.vcf")        
            
            output_directory = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.vcf")
            
            job = bcbio_variation.tumor_pair_ensemble(input_mutect2, input_vardict, input_samtools, output_directory)
            job.name = "bcbio_ensemble." + tumor_pair.name

            jobs.append(job)

        return jobs

    def combine_tumor_pairs_callers(self):
        """
        Merge snvs and indels
        """

        jobs = []

        input_merged_vcfs = [ os.path.join("pairedVariants", tumor_pair.name , tumor_pair.name + ".ensemble.vcf") for tumor_pair in self.tumor_pairs.itervalues()]
        output = os.path.join("pairedVariants", "allPairs.vcf.gz")

        job = concat_jobs([
            gatk.combine_variants(input_merged_vcfs, output)],
            name="combine_tumor_pairs.allPairs")
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
            self.paired_mutect2,
            self.merge_mutect2,
            self.samtools_paired,
            self.merge_filter_paired_samtools,
            self.vardict_paired,
            self.merge_paired_vardict,
            self.ensemble_tumor_pairs,
            self.combine_tumor_pairs_callers
        ]

if __name__ == '__main__': 
    TumorPair()
