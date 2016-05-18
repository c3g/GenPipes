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
from bfx import htslib
from bfx import vcflib
from bfx import vt
from bfx import snpeff
from bfx import gemini
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
                    gatk.mutect2(inputNormal, inputTumor, os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz"))
                ], name="gatk_mutect2." + tumor_pair.name))

            else:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

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
            output = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf.gz")
            output_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf.gz")
            output_somatic =  os.path.join(pair_directory, tumor_pair.name + ".mutect2.somatic.vcf.gz")

            if nb_jobs == 1:
                input = os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz")
                jobs.append(concat_jobs([
                    Job([input], [output], command="ln -s -f " + input + " " + output), 
                    pipe_jobs([
                        Job([output_gz], [None], command="zcat " + output + " | sed 's/TUMOR/"+ tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"+ tumor_pair.normal.name + "/g' "),
                        bcftools.pass_filter(None, None),
                        htslib.bgzip_tabix_vcf(None, output_somatic),
                    ]),
                ], name="symlink_mutect_vcf." + tumor_pair.name))
    
            elif nb_jobs > 1:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                # Create one separate job for each of the first sequences
                inputVCFs = []
                for idx,sequences in enumerate(unique_sequences_per_job):
                    inputVCFs.append(os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".mutect2.vcf.gz"))
                inputVCFs.append(os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"))

                jobs.append(concat_jobs([
                    gatk.cat_variants(inputVCFs, output),
                    pipe_jobs([
                        Job([output_gz], [None], command="zcat " + output + " | sed 's/TUMOR/"+ tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"+ tumor_pair.normal.name + "/g' "),
                        bcftools.pass_filter(None, None),
                        htslib.bgzip_tabix_vcf(None, output_somatic),
                     ]),
                 ], name = "gatk_cat_mutect." + tumor_pair.name))

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
                    ]),
                ], name="samtools_paired." + tumor_pair.name))

            else:
                for region in self.generate_approximate_windows(nb_jobs): #for idx,sequences in enumerate(unique_sequences_per_job):
                    jobs.append(concat_jobs([
                        Job(command="mkdir -p " + samtools_directory),
                        pipe_jobs([
                            samtools.mpileup(paired_sample, None, config.param('samtools_paired', 'mpileup_other_options'), region),
                            samtools.bcftools_view("-", os.path.join(samtools_directory,  tumor_pair.name + "." + region + ".bcf"), config.param('samtools_paired', 'bcftools_view_options'), pair_calling=True),
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
        """

        jobs = []
        nb_jobs = config.param('merge_filter_paired_samtools', 'approximate_nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            samtools_directory = os.path.join(pair_directory, "rawSamtools")
            output = os.path.join(samtools_directory,  tumor_pair.name + ".samtools.bcf")
            output_vcf = os.path.join(pair_directory,  tumor_pair.name + ".samtools.vcf.gz")
            output_somatics = os.path.join(pair_directory,  tumor_pair.name + ".samtools.somatic.vcf.gz")

            if nb_jobs == 1:
                inputs = os.path.join(samtools_directory,  tumor_pair.name + ".bcf") 
                jobs.append(concat_jobs([
                    samtools.bcftools_cat(inputs, output),
                    pipe_jobs([
                        samtools.bcftools_view(output, None),
                        htslib.bgzip_tabix_vcf(None, output_vcf),
                    ]),
                    pipe_jobs([
                        vcflib.vcfsamplediff(tumor_pair.normal.name, tumor_pair.tumor.name, output_vcf, None),
                        bcftools.filter(None, None, config.param('merge_filter_paired_samtools', 'samtools_filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_somatics),
                    ]),
                ], name = "merge_filter_paired_samtools." + tumor_pair.name))

            else:
                inputs = [os.path.join(samtools_directory,  tumor_pair.name + "." + region + ".bcf") for region in self.generate_approximate_windows(nb_jobs)] #for idx,sequences in enumerate(unique_sequences_per_job):
                jobs.append(concat_jobs([
                    samtools.bcftools_cat(inputs, output),
                    pipe_jobs([
                        samtools.bcftools_view(output, None),
                        htslib.bgzip_tabix_vcf(None, output_vcf),
                    ]),
                    pipe_jobs([
                        vcflib.vcfsamplediff(tumor_pair.normal.name, tumor_pair.tumor.name, output_vcf, None),
                        bcftools.filter(None, None, config.param('merge_filter_paired_samtools', 'somatic_filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_somatics),   
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

        jobs = []

        nb_jobs = config.param('vardict_paired', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of vardict jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        use_bed = config.param('vardict_paired', 'use_bed', type='boolean', required=True)
        genome_dictionary = config.param('DEFAULT', 'genome_dictionary', type='filepath')
        if use_bed :
            bed = self.samples[0].readsets[0].beds[0]
            bed_intervals, interval_size = bed_file.parse_bed_file(bed)
            bed_file_list = bed_file.split_by_size(bed_intervals, interval_size, nb_jobs, output="./vardict.tmp")

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            mkdir_job = Job(command="mkdir -p " + vardict_directory)

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
                            Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, 'N', $4) } {print}'"),
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
                            Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, 'N', $4) } {print}'"),
                            htslib.bgzip_tabix_vcf(None, output),
                        ]),
                     ], name="vardict_paired." + tumor_pair.name + ".0"))
                else:
                    bedJob = vardict.dict2beds(genome_dictionary, beds)
                    jobs.append(concat_jobs([mkdir_job,bedJob], name="vardict.genome.beds." + tumor_pair.name))

                    for idx in range(nb_jobs):
                        output = os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz")
                        jobs.append(concat_jobs([
                            mkdir_job,
                            pipe_jobs([
                                vardict.paired_java(inputNormal, inputTumor, tumor_pair.name, None, beds[idx]),
                                vardict.testsomatic(None, None),
                                vardict.var2vcf(None, tumor_pair.normal.name, tumor_pair.tumor.name, None),
                                Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, 'N', $4) } {print}'"),
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
            output = os.path.join(pair_directory,  tumor_pair.name + ".vardict.vcf.gz")
            output_somatic = os.path.join(pair_directory,  tumor_pair.name + ".vardict.somatic.vcf.gz")

            if nb_jobs == 1:
                inputs = os.path.join(vardict_directory,  tumor_pair.name + ".0.vardict.vcf.gz")
                jobs.append(concat_jobs([
                    Job([input], [output], command="ln -s -f " + input + " " + output),
                    pipe_jobs([
                       bcftools.view(output, None, config.param('merge_paired_vardict', 'somatic_filter_options')),
                       htslib.bgzip_tabix_vcf(None, output_somatics),
                    ]),
                ], removable_files=[output], name="symlink_vardict_vcf." + tumor_pair.name ))
            else:
                inputVCFs = []
                for idx in range(nb_jobs):
                    inputVCFs.append(os.path.join(vardict_directory,  tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz"))
                
                jobs.append(concat_jobs([
                    #gatk.cat_variants(inputVCFs, output),
                    pipe_jobs([
                        bcftools.concat(inputVCFs, None),
                        htslib.bgzip_tabix_vcf(None, output),
                    ]),
                    pipe_jobs([
                       bcftools.view(output, None, config.param('merge_filter_paired_vardict', 'somatic_filter_options')),
                       htslib.bgzip_tabix_vcf(None, output_somatic),
                    ]),
                ], name = "merge_filter_paired_vardict." + tumor_pair.name ))

        return jobs


    def ensemble_tumor_pairs(self):
        """
        Apply Bcbio.variations ensemble approach for mutect2, Vardict and Samtools call
        Filter ensemble calls to retain only calls overlapping 2 or more callers
        """

        jobs = []
        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_directory = os.path.join("pairedVariants", tumor_pair.name)
            input_mutect2 = os.path.join(input_directory, tumor_pair.name + ".mutect2.somatic.vcf.gz")
            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.somatic.vcf.gz")
            input_samtools = os.path.join(input_directory, tumor_pair.name + ".samtools.somatic.vcf.gz")        
                    
            output_ensemble = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.vcf")
            output_gz = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.vcf.gz")
            output_flt = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.flt.vcf.gz")

            mkdir_job = Job(command="mkdir -p " + paired_ensemble_directory)   
       
            jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    bcbio_variation.tumor_pair_ensemble(input_mutect2, input_vardict, input_samtools, output_ensemble),
                    htslib.bgzip_tabix_vcf(output_ensemble, output_gz),
                    pipe_jobs([
                        Job([output_gz], [None], command="zgrep -Pv \"set=samtools\\t\" " + output_gz + " | grep -Pv \"set=vardict\\t\" "),
                        htslib.bgzip_tabix_vcf(None, output_flt),
                    ]),
                ], name="bcbio_ensemble." + tumor_pair.name))

        return jobs


    def gatk_variant_annotator(self):
        """
        Add vcf annotations to ensemble vcf: most importantly the AD field
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            inputNormal = os.path.join("alignment",tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment",tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            input_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.flt.vcf.gz")
            output_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.annot.vcf.gz")

            job = gatk.variant_annotator( inputNormal, inputTumor, input_variants, output_variants)
            job.name = "gatk_variant_annotator." + tumor_pair.name            
        
            jobs.append(job)
        
        return jobs

    def compute_cancer_effects(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input = os.path.join(paired_directory, tumor_pair.name + ".ensemble.annot.vcf.gz")
            output = os.path.join(paired_directory, tumor_pair.name + ".ensemble.annot.snpeff.vcf")         
            output_gz = os.path.join(paired_directory, tumor_pair.name + ".ensemble.annot.snpeff.vcf.gz")

            cancer_pair_filename = os.path.join(tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pairs = ("tumor_pair.normal.name\ttumor_pair.tumor.name\n")
            cancer_pair.write(cancer_pairs)

            jobs.append(concat_jobs([           
                snpeff.compute_cancer_effects(input, cancer_pair_filename, output),
                htslib.bgzip_tabix_vcf(output, output_gz)      
            ],name = "compute_cancer_effects." + tumor_pair.name))

        return jobs

    def combine_tumor_pairs_callers(self):
        """
        Combine numerous ensemble vcfs into one vcf for gemini annotations
        """

        jobs = []
        
        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input_merged_vcfs = [os.path.join(os.path.abspath(ensemble_directory), tumor_pair.name , tumor_pair.name + ".ensemble.annot.vcf.gz") for tumor_pair in self.tumor_pairs.itervalues()]
        output = os.path.join(ensemble_directory, "allPairs.ensemble.vcf.gz")

        if len(input_merged_vcfs) == 1:
            job = Job(input_merged_vcfs, [output], command="ln -s -f " + input_merged_vcfs.pop() + " " + output)
            job.name="combine_tumor_pairs.allPairs"
            jobs.append(job)

        else:
            
            jobs.append(concat_jobs([
                gatk.combine_variants(input_merged_vcfs, output)
            ], name="combine_tumor_pairs.allPairs"))

        return jobs

    def decompose_and_normalize_mnps(self):
        """
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt)
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.vt.vcf.gz")

        job = vt.decompose_and_normalize_mnps(input, output)
        job.name = "decompose_and_normalize_mnps.allPairs"

        jobs.append(job)

        return jobs

    def all_pairs_compute_cancer_effect(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        Applied to all tumor pairs.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input =  os.path.join(ensemble_directory, "allPairs.ensemble.vt.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.vt.snpeff.vcf")
        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.vt.snpeff.vcf.gz")

        cancer_pair_filename = os.path.join('cancer_snpeff.tsv')
        cancer_pair = open(cancer_pair_filename, 'w')

        for tumor_pair in self.tumor_pairs.itervalues():
            cancer_pairs = (tumor_pair.normal.name +"\t" + tumor_pair.tumor.name + "\n")
            cancer_pair.write(cancer_pairs)

        jobs.append(concat_jobs([
            snpeff.compute_cancer_effects(input, cancer_pair_filename, output),
            htslib.bgzip_tabix_vcf(output, output_gz),
        ],name = "compute_cancer_effects." + tumor_pair.name))

        return jobs

    def gemini_annotations(self):
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
            Job(command="mkdir -p " + ensemble_directory),
            gemini.gemini_annotations( gemini_prefix + ".ensemble.vt.snpeff.vcf.gz", gemini_prefix + ".gemini." + gemini_version + ".db", temp_dir)
        ], name="gemini_annotations.allPairs"))

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
            self.merge_filter_paired_vardict,
            self.ensemble_tumor_pairs,
            self.gatk_variant_annotator,
            self.compute_cancer_effects,
            self.combine_tumor_pairs_callers,
            self.decompose_and_normalize_mnps,
            self.all_pairs_compute_cancer_effect,
            self.gemini_annotations
        ]

if __name__ == '__main__': 
    TumorPair()
