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
from core.config import config
from core.job import Job, concat_jobs, pipe_jobs

from bfx import bvatools
from bfx import gq_seq_utils
from bfx import gatk
from bfx import igvtools
from bfx import picard
from bfx import samtools
from bfx import tools
from bfx import varscan
from bfx import htslib
from bfx import vt
from bfx import snpeff
from bfx import gemini
from pipelines.dnaseq import dnaseq

log = logging.getLogger(__name__)

class DnaSeqHighCoverage(dnaseq.DnaSeq):
    """
    DnaSeq high Coverage Pipeline
    =================

    The DnaSeq high Coverage Pipeline is based on the DNA-Seq pipeline and follow the first initial steps.
    The difference starts at the Mark Duplicates step. Since this data is high coverage Mark Dup is not run.
    Recalibration is not run either because typically, these datasets are targetted with amplicons or custom
    capture which render recalibration useless.

    Also variant calling is done only using frequency. Not Bayesian callers are used because these typically
    don't fare well with the high coverage.
    """

    def __init__(self, protocol=None):
        self._protocol=protocol
        # Add pipeline specific arguments
        super(DnaSeqHighCoverage, self).__init__(protocol)

    def picard_fixmate(self):
        """
        Verify mate-pair information between mates and fix if needed.
        This ensures that all mate-pair information is in sync between each read and its mate pair.
        Fix is done using [Picard](http://broadinstitute.github.io/picard/).
        """
        jobs = []
        for sample in self.samples:
            sample_directory = os.path.join("alignment", sample.name)
            input_file = os.path.join(sample_directory, sample.name + ".realigned.sorted.bam")
            output_file = os.path.join(sample_directory, sample.name + ".matefixed.sorted.bam")

            job = picard.fix_mate_information(input_file, output_file)
            job.name = "picard_fix_mate_information." + sample.name
            job.samples = [sample]
            jobs.append(job)

        return jobs

    def metrics(self):
        """
        Compute metrics and generate coverage tracks per sample. Multiple metrics are computed at this stage:
        Number of raw reads, Number of filtered reads, Number of aligned reads, Number of duplicate reads,
        Median, mean and standard deviation of insert sizes of reads after alignment, percentage of bases
        covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads)
        whole genome or targeted percentage of bases covered at X reads (%_bases_above_50 means the % of exons
        bases which have at least 50 reads). A TDF (.tdf) coverage track is also generated at this step
        for easy visualization of coverage in the IGV browser.
        """

        jobs = []
        for sample in self.samples:
            input_file_prefix = os.path.join("alignment", sample.name, sample.name + ".matefixed.sorted.")
            input = input_file_prefix + "bam"

            job = picard.collect_multiple_metrics(input, input_file_prefix + "all.metrics")
            job.name = "picard_collect_multiple_metrics." + sample.name
            job.samples = [sample]
            jobs.append(job)

            # Compute genome or target coverage with BVATools
            job = bvatools.depth_of_coverage(
                input,
                input_file_prefix + "coverage.tsv",
                bvatools.resolve_readset_coverage_bed(sample.readsets[0]),
                other_options=config.param('bvatools_depth_of_coverage', 'other_options', required=False)
            )
            job.name = "bvatools_depth_of_coverage." + sample.name
            job.samples = [sample]
            jobs.append(job)

            job = igvtools.compute_tdf(input, input + ".tdf")
            job.name = "igvtools_compute_tdf." + sample.name
            job.samples = [sample]
            jobs.append(job)

        return jobs

    def picard_calculate_hs_metrics(self):
        """
        Compute on target percent of hybridisation based capture.
        """

        jobs = []

        created_interval_lists = []

        for sample in self.samples:
            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            if coverage_bed:
                interval_list = re.sub("\.[^.]+$", ".interval_list", os.path.basename(coverage_bed))

                if not interval_list in created_interval_lists:
                    job = tools.bed2interval_list(None, coverage_bed, interval_list)
                    job.name = "interval_list." + os.path.basename(coverage_bed)
                    jobs.append(job)
                    created_interval_lists.append(interval_list)

                input_file_prefix = os.path.join("alignment", sample.name, sample.name + ".matefixed.sorted.")
                job = picard.calculate_hs_metrics(input_file_prefix + "bam", input_file_prefix + "onTarget.tsv", interval_list)
                job.name = "picard_calculate_hs_metrics." + sample.name
                job.samples = [sample]
                jobs.append(job)
        return jobs

    def gatk_callable_loci(self):
        """
        Computes the callable region or the genome as a bed track.
        """

        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")

            job = gatk.callable_loci(alignment_file_prefix + "matefixed.sorted.bam", alignment_file_prefix + "callable.bed", alignment_file_prefix + "callable.summary.txt")
            job.name = "gatk_callable_loci." + sample.name
            job.samples = [sample]
            jobs.append(job)

        return jobs

    def call_variants(self):
        """
        VarScan caller for insertions and deletions.
        """

        jobs = []

        nb_jobs = config.param('varscan', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of VarScan jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        variants_directory = os.path.join("variants")
        varscan_directory = os.path.join(variants_directory, "rawVarScan")

        beds = []
        for idx in range(nb_jobs):
            beds.append(os.path.join(varscan_directory, 'chrs.' + str(idx) + '.bed'))

        genome_dictionary = config.param('DEFAULT', 'genome_dictionary', type='filepath')

        if nb_jobs > 1:
            bedJob = tools.dict2beds(genome_dictionary, beds)
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + varscan_directory),
                bedJob
            ], name="varscan.genome.beds"))

        bams=[]
        sampleNamesFile = 'varscan_samples.tsv'
        sampleNames = open(sampleNamesFile, 'w')

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            input = os.path.join(alignment_directory, sample.name + ".matefixed.sorted.bam")
            bams.append(input)
            sampleNames.write("%s\n" % sample.name)
            bedfile = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            #sampleNames.append(sample.name)

        if nb_jobs == 1:
            job = concat_jobs([
                Job(command="mkdir -p " + varscan_directory, samples=self.samples),
                pipe_jobs([
                    samtools.mpileup(bams, None, config.param('varscan', 'mpileup_other_options'), regionFile=bedfile),
                    varscan.mpileupcns(None, None, sampleNamesFile, config.param('varscan', 'other_options')),
                    htslib.bgzip_tabix(None, os.path.join(variants_directory, "allSamples.vcf.gz"))
                ])
            ], name="varscan.single")

            jobs.append(job)

        else:
            output_vcfs=[]
            for idx in range(nb_jobs):
                output_vcf = os.path.join(varscan_directory, "allSamples."+str(idx)+".vcf.gz")
                varScanJob = pipe_jobs([
                    samtools.mpileup(bams, None, config.param('varscan', 'mpileup_other_options'), regionFile=beds[idx]),
                    varscan.mpileupcns(None, None, sampleNamesFile, config.param('varscan', 'other_options')),
                    htslib.bgzip_tabix(None, output_vcf)
                ], name = "varscan." + str(idx))
                varScanJob.samples = self.samples
                output_vcfs.append(output_vcf)
                jobs.append(varScanJob)

            job = gatk.cat_variants(output_vcfs, os.path.join(variants_directory, "allSamples.vcf.gz"))
            job.name="gatk_cat_varscan"
            job.samples = self.samples
            jobs.append(job)
        return jobs

    def preprocess_vcf(self):
        """
        Preprocess vcf for loading into a annotation database - gemini : http://gemini.readthedocs.org/en/latest/index.html
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt) and 
        vcf FORMAT modification for correct loading into gemini
        """

        jobs = []

        output_directory = "variants"
        prefix = os.path.join(output_directory, "allSamples")
        outputPreprocess = prefix + ".prep.vt.vcf.gz"
        outputFix = prefix + ".prep.vt.fix.vcf.gz"

        jobs.append(concat_jobs([
            tools.preprocess_varscan( prefix + ".vcf.gz",  prefix + ".prep.vcf.gz" ), 
            pipe_jobs([
                    vt.decompose_and_normalize_mnps(prefix + ".prep.vcf.gz", None),
                    htslib.bgzip_tabix(None, prefix + ".prep.vt.vcf.gz"),
            ]),
            Job(
                [outputPreprocess],
                [outputFix],
                command="zless " + outputPreprocess + " | grep -v 'ID=AD_O' | awk ' BEGIN {OFS=\"\\t\"; FS=\"\\t\"} {if (NF > 8) {for (i=9;i<=NF;i++) {x=split($i,na,\":\") ; if (x > 1) {tmp=na[1] ; for (j=2;j<x;j++){if (na[j] == \"AD_O\") {na[j]=\"AD\"} ; if (na[j] != \".\") {tmp=tmp\":\"na[j]}};$i=tmp}}};print $0} ' | bgzip -cf >  " + outputFix,
                samples=self.samples
            ),
            tools.preprocess_varscan( outputFix,  prefix + ".vt.vcf.gz" )
        ], name="preprocess_vcf.allSamples"))

        return jobs

    def snp_effect(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        """

        jobs = []

        output_directory = "variants"
        snpeff_prefix = os.path.join(output_directory, "allSamples")

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + output_directory, samples=self.samples),
            snpeff.compute_effects( snpeff_prefix + ".vt.vcf.gz", snpeff_prefix + ".vt.snpeff.vcf", split=True),
            htslib.bgzip_tabix( snpeff_prefix + ".vt.snpeff.vcf", snpeff_prefix + ".vt.snpeff.vcf.gz")            
        ], name="compute_effects.allSamples"))

        return jobs

    def gemini_annotations(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        output_directory = "variants"
        temp_dir = os.path.join(os.getcwd(), output_directory)
        gemini_prefix = os.path.join(output_directory, "allSamples")
        gemini_module=config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2],gemini_module[-1]])

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + output_directory, samples=self.samples),
            gemini.gemini_annotations( gemini_prefix + ".vt.snpeff.vcf.gz", gemini_prefix + ".gemini."+gemini_version+".db", temp_dir)
        ], name="gemini_annotations.allSamples"))

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
            self.picard_fixmate,
            self.metrics,
            self.picard_calculate_hs_metrics,
            self.gatk_callable_loci,
            self.call_variants,
            self.preprocess_vcf,
            self.snp_effect,
            self.gemini_annotations,
            self.cram_output
        ]

if __name__ == '__main__': 
    DnaSeqHighCoverage()
