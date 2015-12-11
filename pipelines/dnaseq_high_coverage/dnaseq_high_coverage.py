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
from bfx.sequence_dictionary import *

from bfx import bvatools
from bfx import gq_seq_utils
from bfx import gatk
from bfx import igvtools
from bfx import picard
from bfx import samtools
from bfx import tools
from bfx import varscan
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

    def __init__(self):
        # Add pipeline specific arguments
        super(DnaSeqHighCoverage, self).__init__()

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
            jobs.append(job)

            # Compute genome or target coverage with BVATools
            job = bvatools.depth_of_coverage(
                input,
                input_file_prefix + "coverage.tsv",
                bvatools.resolve_readset_coverage_bed(sample.readsets[0]),
                other_options=config.param('bvatools_depth_of_coverage', 'other_options', required=False)
            )

            job.name = "bvatools_depth_of_coverage." + sample.name
            jobs.append(job)

            job = igvtools.compute_tdf(input, input + ".tdf")
            job.name = "igvtools_compute_tdf." + sample.name
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
                interval_list = re.sub("\.[^.]+$", ".interval_list", coverage_bed)

                if not interval_list in created_interval_lists:
                    job = tools.bed2interval_list(None, coverage_bed, interval_list)
                    job.name = "interval_list." + os.path.basename(coverage_bed)
                    jobs.append(job)
                    created_interval_lists.append(interval_list)

                input_file_prefix = os.path.join("alignment", sample.name, sample.name + ".matefixed.sorted.")
                job = picard.calculate_hs_metrics(input_file_prefix + "bam", input_file_prefix + "onTarget.tsv", interval_list)
                job.name = "picard_calculate_hs_metrics." + sample.name
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
        mkdir_job = Job(command="mkdir -p " + varscan_directory)
        if nb_jobs > 1:
            bedJob = tools.dict2beds(genome_dictionary, beds)
            jobs.append(concat_jobs([mkdir_job,bedJob], name="varscan.genome.beds"))

        bams=[]
        sampleNames=[]
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            input = os.path.join(alignment_directory, sample.name + ".matefixed.sorted.bam")
            bams.append(input)
            sampleNames.append(sample.name)

        if nb_jobs == 1:
            # Create output directory since it is not done by default by GATK tools
            jobs.append(pipe_jobs([
                    samtools.mpitleup(bams, None, config.param('varscan', 'mpileup_other_options')),
                    varscan.mpileupcns(None, os.path.join(variants_directory, "allCalls.vcf"), sampleNames)
                ], name="varscan"))
        else:
            output_vcfs=[]
            for idx in range(nb_jobs):
                output_vcf = os.path.join(varscan_directory, "allCalls."+str(idx)+".vcf")
                varScanJob = pipe_jobs([
                    samtools.mpileup(bams, None, config.param('varscan', 'mpileup_other_options'), regionFile=beds[idx]),
                    varscan.mpileupcns(None, output_vcf, sampleNames)
                ], name = "varscan." + str(idx))
                output_vcfs.append(output_vcf)
                jobs.append(varScanJob)

            job=gatk.cat_variants(output_vcfs, os.path.join(variants_directory, "allCalls.vcf"))
            job.name="gatk_cat_varscan"
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
            self.metrics,
            self.picard_calculate_hs_metrics,
            self.gatk_callable_loci,
            self.call_variants
        ]

if __name__ == '__main__': 
    DnaSeqHighCoverage()
