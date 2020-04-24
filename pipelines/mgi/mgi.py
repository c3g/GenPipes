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
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs

from pipelines.dnaseq import dnaseq

# from bfx import bwa
# from bfx import bvatools
from bfx import cutadapt
from bfx import fgbio
# from bfx import gatk4
# from bfx import igvtools
from bfx import multiqc
# from bfx import sambamba
# from bfx import bash_cmd as bash

log = logging.getLogger(__name__)


class MGISeq(dnaseq.DnaSeqRaw):
    """
    MGI-Seq Pipeline
    ================

    pwet
    """

    def __init__(self, protocol=None):
        self._protocol = protocol
        # Add pipeline specific arguments
        super(MGISeq, self).__init__(protocol)


    def cutadapt(self):
        """
        Raw reads quality trimming and removing of MGI adapters is performed using [Cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html).
        'Adapter1' and 'Adapter2' columns from the readset file ar given to Cutadapt. For PAIRED_END readsets, both adapters are used.
        For SINGLE_END readsets, only Adapter1 is used and left unchanged.
        To trim the front of the read use adapter_5p_fwd and adapter_5p_rev (for PE only) in cutadapt section of ini file.

        This step takes as input files:

        1. FASTQ files from the readset file if available
        2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """

        jobs = []

        for readset in self.readsets:
            output_dir = os.path.join("trim", readset.sample.name)
            trim_file_prefix = os.path.join(output_dir, readset.name)

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[readset.fastq1, readset.fastq2]]
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                adapter_fwd = readset.adapter1
                adapter_rev = readset.adapter2

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[readset.fastq1]]
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                adapter_fwd = readset.adapter1
                adapter_rev = None

            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            jobs.append(
                concat_jobs([
                    Job(
                        command="mkdir -p " + output_dir,
                        removable_files=[output_dir],
                        samples=[readset.sample]
                        ),
                    cutadapt.trim(
                        fastq1,
                        fastq2,
                        trim_file_prefix,
                        adapter_fwd,
                        adapter_rev
                        )
                    ],
                    name="cutadapt." + readset.name)
                )

        return jobs


    # def bwa_mem_sambamba_sort_sam(self):
    #     """
    #     The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
    #     The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
    #     BWA output BAM files are then sorted by coordinate using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

    #     This step takes as input files:

    #     1. Trimmed FASTQ files if available
    #     2. Else, FASTQ files from the readset file if available
    #     3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
    #     """

    #     jobs = []
    #     for readset in self.readsets:
    #         trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
    #         alignment_directory = os.path.join("alignment", readset.sample.name)
    #         readset_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam")
    #         index_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam.bai")

    #         # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
    #         if readset.run_type == "PAIRED_END":
    #             candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
    #             if readset.fastq1 and readset.fastq2:
    #                 candidate_input_files.append([readset.fastq1, readset.fastq2])
    #             if readset.bam:
    #                 prefix = os.path.join(
    #                     self.output_dir,
    #                     "raw_reads",
    #                     readset.sample.name,
    #                     re.sub("\.bam$", ".", os.path.basename(readset.bam))
    #                 )
    #                 candidate_input_files.append([prefix + "pair1.fastq.gz", prefix + "pair2.fastq.gz"])
    #             [fastq1, fastq2] = self.select_input_files(candidate_input_files)

    #         elif readset.run_type == "SINGLE_END":
    #             candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
    #             if readset.fastq1:
    #                 candidate_input_files.append([readset.fastq1])
    #             if readset.bam:
    #                 prefix = os.path.join(
    #                     self.output_dir,
    #                     "raw_reads",
    #                     readset.sample.name,
    #                     re.sub("\.bam$", ".", os.path.basename(readset.bam))
    #                 )
    #                 candidate_input_files.append([prefix + ".single.fastq.gz"])
    #             [fastq1] = self.select_input_files(candidate_input_files)
    #             fastq2 = None

    #         else:
    #             _raise(SanitycheckError("Error: run type \"" + readset.run_type +
    #             "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

    #         jobs.append(
    #             concat_jobs([
    #                 bash.mkdir(os.path.dirname(readset_bam)),
    #                 pipe_jobs([
    #                     bwa.mem(
    #                         fastq1,
    #                         fastq2,
    #                         read_group="'@RG" + \
    #                             "\\tID:" + readset.name + \
    #                             "\\tSM:" + readset.sample.name + \
    #                             "\\tLB:" + (readset.library if readset.library else readset.sample.name) + \
    #                             ("\\tPU:run" + readset.run + "_" + readset.lane if readset.run and readset.lane else "") + \
    #                             ("\\tCN:" + config.param('bwa_mem', 'sequencing_center') if config.param('bwa_mem', 'sequencing_center', required=False) else "") + \
    #                             "\\tPL:MGI" + \
    #                             "'"
    #                         ),
    #                     sambamba.view(
    #                         "/dev/stdin",
    #                         None,
    #                         "-S -f bam"
    #                         ),
    #                     sambamba.sort(
    #                         "/dev/stdin",
    #                         readset_bam,
    #                         config.param('sambamba_sort_sam', 'tmp_dir', required=True)
    #                         )
    #                     ]),
    #                 sambamba.index(
    #                     readset_bam,
    #                     index_bam
    #                     )
    #                 ],
    #                 name="bwa_mem_sambamba_sort_sam." + readset.name,
    #                 samples=[readset.sample]
    #                 )
    #             )

    #     return jobs

    def fgbio_trim_primers(self):
        """
        Remove primer sequences to individual bam files using fgbio
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            input_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            output_bam = os.path.join(alignment_directory, sample.name + ".sorted.primerTrim.bam")

            jobs.append(
                concat_jobs([
                    Job(command="mkdir -p " + os.path.dirname(output_bam)),
                    fgbio.trim_primers(
                        input_bam,
                        output_bam,
                        hard_clip=True
                        )
                    ],
                    name="fgbio_trim_primers." + sample.name,
                    samples=[sample]
                    )
                )

        return jobs


    def mgi_metrics(self):
        """
        Multiple metrcis from dnaseq
        """

        jobs = []

        jobs.extend(self.metrics_dna_picard_metrics())
        jobs.extend(self.metrics_dna_sample_qualimap())
        jobs.extend(self.metrics_dna_sambamba_flagstat())
        jobs.extend(self.picard_calculate_hs_metrics())
        jobs.extend(self.metrics())

        return jobs


    def mgi_calling(self):
        """
        Calling steps from dnaseq
        """

        jobs = []

        jobs.extend(self.gatk_haplotype_caller())
        jobs.extend(self.merge_and_call_individual_gvcf())
        # jobs.extend(self.combine_gvcf())
        # jobs.extend(self.merge_and_call_combined_gvcf())
        # jobs.extend(self.variant_recalibrator())
        # jobs.extend(self.haplotype_caller_decompose_and_normalize())
        # jobs.extend(self.haplotype_caller_flag_mappability())
        # jobs.extend(self.haplotype_caller_snp_id_annotation())
        # jobs.extend(self.haplotype_caller_snp_effect())
        # jobs.extend(self.haplotype_caller_dbnsfp_annotation())
        # jobs.extend(self.haplotype_caller_gemini_annotations())
        # jobs.extend(self.haplotype_caller_metrics_vcf_stats())

        return jobs

    def run_multiqc(self):

        jobs = []

        metrics_directory = os.path.join("metrics", "dna")
        input_dep = []
        inputs = []
        for sample in self.samples:
            input_oxog = os.path.join(metrics_directory, sample.name, "picard_metrics", sample.name + ".oxog_metrics.txt")
            input_qcbias = os.path.join(metrics_directory, sample.name, "picard_metrics", sample.name + ".qcbias_metrics.txt")
            input_all_picard = os.path.join(metrics_directory, sample.name, "picard_metrics", sample.name + ".all.metrics.quality_distribution.pdf")
            # input_qualimap = os.path.join(metrics_directory, sample.name, "qualimap", sample.name + ".genome_results.txt")
            # input_fastqc = os.path.join(metrics_directory, sample.name, "fastqc", sample.name + ".fastqc.zip")
            # input_flagstat = os.path.join(metrics_directory, sample.name, "flagstat", sample.name + ".flagstat")

            input_dep += [
                input_oxog,
                input_qcbias,
                input_all_picard
                # input_qualimap,
                # input_fastqc,
                # input_flagstat
            ]

            inputs += [os.path.join(metrics_directory, sample.name)]

        output = os.path.join(metrics_directory, "multiqc_report")

        job = multiqc.run(
            inputs,
            output,
            input_dep
            )
        job.name = "multiqc_all_samples"
        job.samples = self.samples

        jobs.append(job)

        return jobs


    @property
    def steps(self):
        return [
            self.cutadapt,
            self.bwa_mem_sambamba_sort_sam,
            self.sambamba_merge_sam_files,
            self.fgbio_trim_primers,
            self.mgi_metrics,
            self.mgi_calling
            # self.run_multiqc
        ]

if __name__ == '__main__':
    MGISeq()
