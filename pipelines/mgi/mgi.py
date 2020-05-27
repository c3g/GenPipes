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
import utils.utils

from pipelines.dnaseq import dnaseq

from bfx import bwa
from bfx import bedtools
from bfx import cutadapt
from bfx import fgbio
# from bfx import gatk4
from bfx import htslib
from bfx import ivar
from bfx import multiqc
from bfx import sambamba
from bfx import samtools
from bfx import bash_cmd as bash

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
            trim_directory = os.path.join("trim", readset.sample.name)
            trim_file_prefix = os.path.join(trim_directory, readset.name)

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
                    bash.mkdir(trim_directory),
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


    def mapping_bwa_mem_sambamba(self):
        """
        The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
        The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
        BWA output BAM files are then sorted by coordinate using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

        This step takes as input files:

        1. Trimmed FASTQ files if available
        2. Else, FASTQ files from the readset file if available
        3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """

        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_directory = os.path.join("alignment", readset.sample.name)
            readset_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam")
            index_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam.bai")

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    prefix = os.path.join(
                        self.trim_directory,
                        "raw_reads",
                        readset.sample.name,
                        re.sub("\.bam$", ".", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([prefix + "pair1.fastq.gz", prefix + "pair2.fastq.gz"])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    prefix = os.path.join(
                        self.trim_directory,
                        "raw_reads",
                        readset.sample.name,
                        re.sub("\.bam$", ".", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([prefix + ".single.fastq.gz"])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None

            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(readset_bam)),
                    pipe_jobs([
                        bwa.mem(
                            fastq1,
                            fastq2,
                            read_group="'@RG" + \
                                "\\tID:" + readset.name + \
                                "\\tSM:" + readset.sample.name + \
                                "\\tLB:" + (readset.library if readset.library else readset.sample.name) + \
                                ("\\tPU:run" + readset.run + "_" + readset.lane if readset.run and readset.lane else "") + \
                                ("\\tCN:" + config.param('mapping_bwa_mem_sambamba', 'sequencing_center') if config.param('mapping_bwa_mem_sambamba', 'sequencing_center', required=False) else "") + \
                                ("\\tPL:" + config.param('mapping_bwa_mem_sambamba', 'sequencing_technology') if config.param('mapping_bwa_mem_sambamba', 'sequencing_technology', required=False) else "Illumina") + \
                                "'",
                                ini_section='mapping_bwa_mem_sambamba',
                                other_options=config.param('mapping_bwa_mem_sambamba', 'bwa_other_options')
                                ),
                        sambamba.view(
                            "/dev/stdin",
                            None,
                            options=config.param('mapping_bwa_mem_sambamba', 'sambamba_view_other_options')
                            ),
                        sambamba.sort(
                            "/dev/stdin",
                            readset_bam,
                            tmp_dir=config.param('mapping_bwa_mem_sambamba', 'tmp_dir', required=True),
                            other_options=config.param('mapping_bwa_mem_sambamba', 'sambamba_sort_other_options', required=False)
                            )
                        ]),
                    sambamba.index(
                        readset_bam,
                        index_bam,
                        other_options=config.param('mapping_bwa_mem_sambamba', 'sambamba_index_other_options', required=False)
                        )
                    ],
                    name="mapping_bwa_mem_sambamba." + readset.name,
                    samples=[readset.sample]
                    )
                )

        return jobs

    def sambamba_filtering(self):
        """
        Filter raw bams with sambamba and an awk cmd to filter by insert size
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            input_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            output_bam = os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(output_bam)),
                    pipe_jobs([
                        sambamba.view(
                            input_bam=input_bam,
                            output_bam=None,
                            options="-h"
                            ),
                        Job(
                            input_files=[],
                            output_files=[],
                            command="""awk 'substr($0,1,1)=="@" ||\\
 ($9 >= {min_insert_size} && $9 <= {max_insert_size}) ||\\
 ($9 <= -{min_insert_size} && $9 >= -{max_insert_size})'""".format(
     min_insert_size=config.param('sambamba_filtering', 'min_insert_size', required=False),
     max_insert_size=config.param('sambamba_filtering', 'max_insert_size', required=False)
     )
                            ),
                        sambamba.view(
                            input_bam="/dev/stdin",
                            output_bam=output_bam,
                            options=config.param('sambamba_filtering', 'sambamba_filtering_other_options', required=False)
                            )
                        ]),
                    ],
                    name="sambamba_filtering." + sample.name,
                    samples=[sample]
                    )
                )

        return jobs


    def fgbio_trim_primers(self):
        """
        Remove primer sequences to individual bam files using fgbio
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            input_bam = os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")
            output_bam = os.path.join(alignment_directory, sample.name + ".sorted.primerTrim.bam")

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(output_bam)),
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


    def ivar_trim_primers(self):
        """
        Remove primer sequences to individual bam files using ivar
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            input_bam = os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")
            output_prefix = os.path.join(alignment_directory, re.sub("\.sorted.filtered.bam$", ".primerTrim", os.path.basename(input_bam)))
            output_bam = os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.bam")

            jobs.append(
                concat_jobs([
                    bash.mkdir(alignment_directory),
                    ivar.trim_primers(
                        input_bam,
                        output_prefix
                        ),
                    sambamba.sort(
                        output_prefix + ".bam",
                        output_bam,
                        config.param('ivar_trim_primers', 'tmp_dir')
                        )
                    ],
                    name="ivar_trim_primers." + sample.name,
                    samples=[sample],
                    removable_files=[output_prefix + ".bam"]
                    )
                )

        return jobs


    def metrics_sambamba_flagstat(self):
        """
        Sambamba flagstsat
        """

        jobs = []
        for sample in self.samples:
            flagstat_directory = os.path.join("metrics", "dna", sample.name, "flagstat")
            alignment_directory = os.path.join("alignment", sample.name)
            [input_bam] = self.select_input_files([
                # [os.path.join(alignment_directory, sample.name + ".sorted.primerTrim.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")],
                # [os.path.join(alignment_directory, sample.name + ".sorted.bam")]
            ])

            # for bam in [os.path.join(alignment_directory, sample.name + ".sorted.primerTrim.bam"), os.path.join(alignment_directory, sample.name + ".sorted.bam")]:
            output = os.path.join(flagstat_directory, re.sub("\.bam$", ".flagstat", os.path.basename(input_bam)))

            jobs.append(
                concat_jobs([
                    bash.mkdir(flagstat_directory),
                    sambamba.flagstat(
                        input_bam,
                        output,
                        config.param('sambamba_flagstat', 'flagstat_options')
                        )
                    ],
                    name="sambamba_flagstat." + sample.name,
                    samples=[sample]
                    )
                )

        return jobs

    def metrics_bedtools_genomecov(self):
        """
        bedtools genome coverage
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            [input_bam] = self.select_input_files([
                # [os.path.join(alignment_directory, sample.name + ".sorted.primerTrim.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.bam")]
            ])

            # for bam in [os.path.join(alignment_directory, sample.name + ".sorted.primerTrim.bam"), os.path.join(alignment_directory, sample.name + ".sorted.bam")]:
            output = os.path.join(alignment_directory, re.sub("\.bam$", ".BedGraph", os.path.basename(input_bam)))

            jobs.append(
                concat_jobs([
                    bash.mkdir(alignment_directory),
                    bedtools.genomecov(
                        input_bam,
                        output
                        )
                    ],
                    name="bedtools_genomecov." + sample.name,#re.sub("\.bam$", "", os.path.basename(bam)),
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
        jobs.extend(self.metrics_sambamba_flagstat())
        jobs.extend(self.metrics_bedtools_genomecov())
        # jobs.extend(self.metrics_dna_sambamba_flagstat())
        jobs.extend(self.picard_calculate_hs_metrics())
        jobs.extend(self.metrics())

        return jobs


    def mgi_calling(self):
        """
        Calling steps from dnaseq
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            # haplotype_directory = os.path.join(alignment_directory, "rawHaplotypeCaller")

            [input_bam] = self.select_input_files([
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.bam")]
            ])

            # for bam in [os.path.join(alignment_directory, sample.name + ".sorted.primerTrim.bam"), os.path.join(alignment_directory, sample.name + ".sorted.bam")]:

            output_prefix = os.path.join(alignment_directory, re.sub("\.bam$", "", os.path.basename(input_bam)))
            output_tsv = output_prefix + ".tsv"
            output_vcf = output_prefix + ".vcf"

            jobs.append(
                concat_jobs([
                    bash.mkdir(alignment_directory),
                    pipe_jobs([
                        samtools.mpileup(
                            input_bams=input_bam,
                            output=None,
                            other_options=config.param('ivar_call_variants', 'mpileup_options'),
                            region=None,
                            regionFile=None,
                            ini_section='ivar_call_variants'
                            ),
                        ivar.call_variants(
                            output_prefix
                            )
                        ]),
                    ivar.tsv_to_vcf(
                        output_tsv,
                        output_vcf
                        ),
                    htslib.bgzip_tabix(
                        output_vcf,
                        output_vcf + ".gz"
                        )
                    ],
                    name="ivar_call_variants." + sample.name,#re.sub("\.bam$", "", os.path.basename(bam)),
                    samples=[sample],
                    removable_files=[output_tsv, output_vcf]
                    )
                )

        # jobs.extend(self.gatk_haplotype_caller())
        # jobs.extend(self.merge_and_call_individual_gvcf())

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

    def ivar_create_consensus(self):
        """
        Remove primer sequences to individual bam files using fgbio
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            [input_bam] = self.select_input_files([
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.bam")]
            ])

            # for bam in [os.path.join(alignment_directory, sample.name + ".sorted.primerTrim.bam"), os.path.join(alignment_directory, sample.name + ".sorted.bam")]:

            output_prefix = os.path.join(alignment_directory, re.sub("\.bam$", ".consensus", os.path.basename(input_bam)))

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(output_prefix)),
                    pipe_jobs([
                        samtools.mpileup(
                            input_bams=input_bam,
                            output=None,
                            other_options=config.param('ivar_create_consensus', 'mpileup_options'),
                            region=None,
                            regionFile=None,
                            ini_section='ivar_create_consensus'
                            ),
                        ivar.create_consensus(
                            output_prefix
                            )
                        ])
                    ],
                    name="ivar_create_consensus." + sample.name,#re.sub("\.bam$", "", os.path.basename(bam)),
                    samples=[sample]
                    )
                )

        return jobs


    def rename_consensus_header(self):
        """
        Rename reads headers
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            [input_fa] = self.select_input_files([
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.consensus.fa")],
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.consensus.fa")],
                [os.path.join(alignment_directory, sample.name + ".sorted.consensus.fa")]
            ])

            output_fa = os.path.join(alignment_directory, re.sub("\.fa$", ".gisaid_renamed.fa", os.path.basename(input_fa)))

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(output_fa)),
                    Job(
                        input_files=[input_fa],
                        output_files=[output_fa],
                        command="""\\
awk '/^>/{{print ">{country}/{province}-{sample}/{year} seq_method:{seq_method}|assemb_method:{assemb_method}|snv_call_method:{snv_call_method}"; next}}{{print}}' < {input_fa} > {output_fa}""".format(
    country=config.param('rename_consensus_header', 'country', required=False),
    province=config.param('rename_consensus_header', 'province', required=False),
    year=config.param('rename_consensus_header', 'year', required=False),
    seq_method=config.param('rename_consensus_header', 'seq_method', required=False),
    assemb_method=config.param('rename_consensus_header', 'assemb_method', required=False),
    snv_call_method=config.param('rename_consensus_header', 'snv_call_method', required=False),
    sample=sample.name,
    input_fa=input_fa,
    output_fa=output_fa
    )
                        )
                ],
                name="rename_consensus_header." + sample.name,
                samples=[sample]
                )
            )

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
            self.mapping_bwa_mem_sambamba,
            self.sambamba_merge_sam_files,
            self.sambamba_filtering,
            # self.fgbio_trim_primers,
            self.ivar_trim_primers,
            self.mgi_metrics,
            self.mgi_calling,
            self.ivar_create_consensus,
            self.rename_consensus_header
            # self.run_multiqc
        ]

if __name__ == '__main__':
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        MGISeq()
