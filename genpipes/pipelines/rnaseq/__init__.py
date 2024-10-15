################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import argparse
import collections
import csv
import logging
import math
import os
import re

# GenPipes Modules
from ...core.config import global_conf, _raise, SanitycheckError
from ...core.job import Job, concat_jobs, pipe_jobs
from .. import common
from ...bfx.sequence_dictionary import parse_sequence_dictionary_file, split_by_size

from ...bfx import (
    adapters,
    annoFuse,
    arriba,
    ballgown,
    bash_cmd as bash,
    bcftools,
    bvatools,
    bwa,
    cpsr,
    deeptools,
    deliverables,
    differential_expression,
    gatk4,
    gemini,
    gtex_pipeline,
    htseq,
    htslib,
    job2json_project_tracking,
    metrics,
    multiqc,
    pcgr,
    picard2 as picard,
    rseqc,
    rnaseqc2,
    sambamba,
    skewer,
    snpeff,
    star,
    star_fusion,
    stringtie,
    tools,
    vcfanno,
    vt,
    sortmerna
    )

log = logging.getLogger(__name__)

class RnaSeqRaw(common.Illumina):
    """
RNA-Seq Pipeline
================

The standard MUGQIC RNA-Seq pipeline has three protocols (stringtie, variants, cancer), stringtie is the
default protocol and applicable in most cases.

All three protocols are based on the use of the [STAR aligner](https://code.google.com/p/rna-star/)
to align reads to the reference genome. These alignments are used during
downstream analysis (for example in stringtie protocol, to determine differential expression of genes and transcripts).

The [StringTie](https://ccb.jhu.edu/software/stringtie/) suite is used for differential transcript expression (DTE) analysis, whereas
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and
[edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) are used for the differential gene expression (DGE) analysis.

The "stringtie" protocol requires a design file which will be used to define the comparison groups
in the differential analyses. The design file format is described [here](https://genpipes.readthedocs.io/en/latest/get-started/concepts/design_file.html).
In addition, [Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html) is used to calculate differential
transcript and gene expression levels and test them for significant differences.
It can also take a batch file (optional) which will be used to correct for batch effects
in the differential analyses. The batch file format is described [here](https://bitbucket.org/mugqic/mugqic_pipelines/src#markdown-header-batch-file)

The variants protocol is used when variant detection, is the main goal of the analysis. GATK best practices workflow
is used for variant calling. It also contains a step for annotating genes using [gemini](https://gemini.readthedocs.io/en/latest/)

The cancer protocol contains all the steps in the variants protocol but it is
specifically designed for cancer data sets due to the
complexity of cancer samples and additional analyses those projects often entail.
The goal of the cancer protocol is comparing expression to known benchmarks. In addition to the steps in the variants
protocol, it contains four specific steps. Three of them (run_star_fusion, run_arriba, run_annofuse)
are related to detection and annotation of gene fusions. For that, [Star-fusion](https://github.com/STAR-Fusion/STAR-Fusion-Tutorial/wiki),
[Arriba](https://arriba.readthedocs.io/en/latest/) and [anno-Fuse](https://rdrr.io/github/d3b-center/annoFuse/) are
used. Furthermore, [RSeQC](http://rseqc.sourceforge.net/) provides RNA-seq quality control metrics to asses the
quality of the data.
    """

    def __init__(self, *args, **kwargs):
        # Add pipeline specific arguments
        super(RnaSeqRaw, self).__init__(*args, **kwargs)

    @property
    def output_dirs(self):
        dirs = {
            'raw_reads_directory': os.path.relpath(os.path.join(self.output_dir, 'raw_reads'), self.output_dir),
            'trim_directory': os.path.relpath(os.path.join(self.output_dir, 'trim'), self.output_dir),            
            'alignment_1stPass_directory': os.path.relpath(os.path.join(self.output_dir, 'alignment_1stPass'), self.output_dir),
            'alignment_directory': os.path.relpath(os.path.join(self.output_dir, 'alignment'), self.output_dir),
            'stringtie_directory': os.path.relpath(os.path.join(self.output_dir, 'stringtie'), self.output_dir),
            'ballgown_directory': os.path.relpath(os.path.join(self.output_dir, 'ballgown'), self.output_dir),
            'DGE_directory': os.path.relpath(os.path.join(self.output_dir, 'DGE'), self.output_dir),
            'raw_counts_directory': os.path.relpath(os.path.join(self.output_dir, 'raw_counts'), self.output_dir),
            'tracks_directory': os.path.relpath(os.path.join(self.output_dir, 'tracks'), self.output_dir),
            'exploratory_directory': os.path.relpath(os.path.join(self.output_dir, 'exploratory'), self.output_dir),
            'metrics_directory': os.path.relpath(os.path.join(self.output_dir, 'metrics'), self.output_dir),
            'report_directory': os.path.relpath(os.path.join(self.output_dir, 'report'), self.output_dir),
            'fusion_directory': os.path.relpath(os.path.join(self.output_dir, 'fusion'), self.output_dir)
        }
        return dirs

    @property
    def multiqc_inputs(self):
        if not hasattr(self, "_multiqc_inputs"):
            self._multiqc_inputs = []
        return self._multiqc_inputs

    @multiqc_inputs.setter
    def multiqc_inputs(self, value):
        self._multiqc_inputs = value

    @property
    def sequence_dictionary(self):
        if not hasattr(self, "_sequence_dictionary"):
            self._sequence_dictionary = parse_sequence_dictionary_file(global_conf.global_get('DEFAULT', 'genome_dictionary', param_type='filepath'), variant=False)
        return self._sequence_dictionary

    def star_genome_length(self):
        """
        Calculation for setting genomeSAindexNbases for STAR index.
        """
        genome_index = csv.reader(open(global_conf.global_get('DEFAULT', 'genome_fasta', param_type='filepath') + ".fai", 'r'), delimiter='\t')

        return int(min(14, math.log2(sum([int(chromosome[1]) for chromosome in genome_index])) / 2 - 1))

    def skewer_trimming(self):
        """
        [Skewer](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-182) is used mainly for
        detection and trimming adapter sequences from raw fastq files. Other features of Skewer is listed
        [here](https://github.com/relipmoc/skewer).
        """

        jobs = []

        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")
        for readset in self.readsets:
            output_dir = os.path.join(self.output_dirs["trim_directory"], readset.sample.name)
            trim_file_prefix = os.path.join(output_dir, readset.name)

            adapter_file = global_conf.global_get('skewer_trimming', 'adapter_file', required=False, param_type='filepath')
            adapter_job = None

            quality_offset = readset.quality_offset

            if not adapter_file:
                adapter_file = os.path.join(output_dir, "adapter.tsv")
                adapter_job = adapters.create(
                    readset,
                    adapter_file
                )

            fastq1 = ""
            fastq2 = ""
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [
                    [
                        readset.fastq1,
                        readset.fastq2
                    ]
                ]
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dirs['raw_reads_directory'],
                        readset.sample.name,
                        re.sub(r"\.bam$", ".", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append(
                        [
                            prefix + "pair1.fastq.gz",
                            prefix + "pair2.fastq.gz"
                        ]
                    )
                    prefix = os.path.join(
                        self.output_dirs['raw_reads_directory'],
                        readset.sample.name,
                        readset.name + "."
                    )
                    candidate_input_files.append(
                        [
                            prefix + "pair1.fastq.gz",
                            prefix + "pair2.fastq.gz"
                        ]
                    )
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[readset.fastq1]]
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dirs['raw_reads_directory'],
                        readset.sample.name,
                        re.sub(r"\.bam$", ".", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([prefix + ".single.fastq.gz"])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None

            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            output_dir,
                            remove=True
                        ),
                        bash.mkdir(
                            link_directory
                            ),
                        adapter_job,
                        skewer.trim(
                            fastq1,
                            fastq2,
                            trim_file_prefix,
                            adapter_file,
                            quality_offset
                        ),
                        bash.ln(
                            os.path.relpath(trim_file_prefix + "-trimmed-pair1.fastq.gz", os.path.dirname(trim_file_prefix + ".trim.pair1.fastq.gz")),
                            trim_file_prefix + ".trim.pair1.fastq.gz",
                            input=trim_file_prefix + "-trimmed-pair1.fastq.gz"
                        ),
                        bash.ln(
                            os.path.relpath(trim_file_prefix + "-trimmed-pair2.fastq.gz", os.path.dirname(trim_file_prefix + ".trim.pair2.fastq.gz")),
                            trim_file_prefix + ".trim.pair2.fastq.gz",
                            input=trim_file_prefix + "-trimmed-pair2.fastq.gz"
                        ),
                        bash.ln(
                            os.path.relpath(trim_file_prefix + "-trimmed.log", link_directory),
                            os.path.join(link_directory, readset.name + "-trimmed.log"),
                            trim_file_prefix + "-trimmed.log"
                            )
                    ],
                    name="skewer_trimming." + readset.name,
                    removable_files=[output_dir],
                    samples=[readset.sample],
                    readsets=[readset]
                )
            )

            self.multiqc_inputs.append(trim_file_prefix + "-trimmed.log")
        return jobs

    def sortmerna(self):
        """
        Calculation of ribosomal RNA per read based on known ribosomal sequences from archea, bacteria and eukaryotes.
        Using [sortmeRNA](https://github.com/sortmerna/sortmerna)

        Taking trimmed fastqs and reporting on each read, either paired-end or single end.
        """

        jobs = []
        for readset in self.readsets:
            output_dir = os.path.join(self.output_dirs["metrics_directory"], "sortmerna")
            output_dir_sample = os.path.join(output_dir, readset.sample.name, readset.name)
            index_directory = os.path.join(output_dir, "idx-dir")
            link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

            trim_fastq1 = ""
            trim_fastq2 = ""

            if readset.run_type == "PAIRED_END":
                trim_fastq1 = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim." + "pair1.fastq.gz")
                trim_fastq2 = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim." + "pair2.fastq.gz")

                sortmerna_job=sortmerna.paired(
                            trim_fastq1,
                            trim_fastq2,
                            output_dir,
                            output_dir_sample,
                            readset.name
                        )
                inputs = [trim_fastq1, trim_fastq2]

            elif readset.run_type == "SINGLE_END":
                trim_fastq1 = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim." + "single.fastq.gz")
                trim_fastq2 = None

                sortmerna_job=sortmerna.single(
                            trim_fastq1,
                            output_dir,
                            output_dir_sample,
                            readset.name
                        )
                inputs = [trim_fastq1]

            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            jobs.append(
                concat_jobs(
                [
                    bash.rm(output_dir_sample),
                    bash.mkdir(output_dir_sample),
                    bash.mkdir(index_directory),
                    bash.mkdir(link_directory),
                    sortmerna_job,
                    bash.ln(
                        os.path.relpath(os.path.join(output_dir_sample, readset.name + ".aligned.log"), link_directory),
                        os.path.join(link_directory, readset.name + ".aligned.log"),
                        os.path.join(output_dir_sample, readset.name + ".aligned.log")
                    )
                ],
                name="sortmerna." + readset.name,
                input_dependency = inputs,
                samples=[readset.sample]
                )
            )

            self.multiqc_inputs.append(os.path.join(output_dir_sample, readset.name + ".aligned.log"))

        return jobs


    def star(self):
        """
        The filtered reads are aligned to a reference genome. The alignment is done per readset of sequencing
        using the [STAR](https://code.google.com/p/rna-star/) software. It generates a Binary Alignment Map file (.bam).

        This step takes as input files:

        1. Trimmed FASTQ files if available
        2. Else, FASTQ files from the readset file if available
        3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """

        jobs = []
        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")
        project_index_directory = "reference.Merged"
        project_junction_file = os.path.join(self.output_dirs["alignment_1stPass_directory"], "AllSamples.SJ.out.tab")
        individual_junction_list=[]
        genome_length = self.star_genome_length()
        mapping = global_conf.global_get("star_align", "mapping", required=False) # option to skip 2-pass mapping with star aligner and only do one pass

        if not mapping or mapping == "2-pass": # if mapping is not explicitely set in config file or is set to 2-pass, start 1st pass alignment, otherwise skip to next pass
            #pass 1 -alignment
            for readset in self.readsets:
                trim_file_prefix = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim.")
                alignment_1stPass_directory = os.path.join(self.output_dirs["alignment_1stPass_directory"], readset.sample.name, readset.name)
                individual_junction_list.append(os.path.join(alignment_1stPass_directory,"SJ.out.tab"))

                if readset.run_type == "PAIRED_END":
                    candidate_input_files = [
                        [
                            trim_file_prefix + "pair1.fastq.gz",
                            trim_file_prefix + "pair2.fastq.gz"
                        ]
                    ]
                    if readset.fastq1 and readset.fastq2:
                        candidate_input_files.append(
                            [
                                readset.fastq1,
                                readset.fastq2
                            ]
                        )
                    if readset.bam:
                        candidate_input_files.append(
                            [
                                re.sub(r"\.bam$", ".pair1.fastq.gz", readset.bam),
                                re.sub(r"\.bam$", ".pair2.fastq.gz", readset.bam)
                            ]
                        )
                    [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                elif readset.run_type == "SINGLE_END":
                    candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                    if readset.fastq1:
                        candidate_input_files.append([readset.fastq1])
                    if readset.bam:
                        candidate_input_files.append([re.sub(r"\.bam$", ".single.fastq.gz", readset.bam)])
                    [fastq1] = self.select_input_files(candidate_input_files)
                    fastq2 = None
                else:
                    _raise(SanitycheckError("Error: run type \"" + readset.run_type + "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

                rg_platform = global_conf.global_get('star_align', 'platform', required=False)
                rg_center = global_conf.global_get('star_align', 'sequencing_center', required=False)

                job = star.align(
                    reads1=fastq1,
                    reads2=fastq2,
                    output_directory=alignment_1stPass_directory,
                    genome_index_folder=None,
                    rg_id=readset.name,
                    rg_sample=readset.sample.name,
                    rg_library=readset.library if readset.library else "",
                    rg_platform_unit=readset.run + "_" + readset.lane if readset.run and readset.lane else "",
                    rg_platform=rg_platform if rg_platform else "",
                    rg_center=rg_center if rg_center else ""
                )
                job.name = "star_align.1." + readset.name
                job.samples = [readset.sample]
                job.readsets = [readset]
                jobs.append(job)

            jobs.append(
                concat_jobs(
                    [
                        #pass 1 - contatenate junction
                        Job(
                            samples=self.samples,
                            readsets=self.readsets
                            ),
                        star.concatenate_junction(
                            input_junction_files_list=individual_junction_list,
                            output_junction_file=project_junction_file
                        ),
                        #pass 1 - genome indexing
                        star.index(
                            genome_index_folder=project_index_directory,
                            junction_file=project_junction_file,
                            genome_length=genome_length
                        )
                    ],
                    name="star_index.AllSamples",
                    samples=self.samples,
                    readsets=self.readsets
                )
            )

        ######
        #Single Pass or Pass 2 - alignment
        #This is the only pass if 1-pass option is specified in config file, pass 2 if not specified or 2-pass option specified
        for readset in self.readsets:
            trim_file_prefix = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim.")
            alignment_pass_directory = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name, readset.name)
            sample_log = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name, readset.name, "Log.final.out")

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [
                    [
                        trim_file_prefix + "pair1.fastq.gz",
                        trim_file_prefix + "pair2.fastq.gz"
                    ]
                ]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append(
                        [
                            readset.fastq1,
                            readset.fastq2
                        ]
                    )
                if readset.bam:
                    candidate_input_files.append(
                        [
                            re.sub(r"\.bam$", ".pair1.fastq.gz", readset.bam),
                            re.sub(r"\.bam$", ".pair2.fastq.gz", readset.bam)
                        ]
                    )
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub(r"\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            rg_platform = global_conf.global_get('star_align', 'platform', required=False)
            rg_center = global_conf.global_get('star_align', 'sequencing_center', required=False)

            job = concat_jobs(
                [
                    star.align(
                    reads1=fastq1,
                    reads2=fastq2,
                    output_directory=alignment_pass_directory,
                    genome_index_folder=project_index_directory if not mapping or mapping=="2-pass" else None,
                    rg_id=readset.name,
                    rg_sample=readset.sample.name,
                    rg_library=readset.library if readset.library else "",
                    rg_platform_unit=readset.run + "_" + readset.lane if readset.run and readset.lane else "",
                    rg_platform=rg_platform if rg_platform else "",
                    rg_center=rg_center if rg_center else "",
                    create_wiggle_track=True,
                    search_chimeres=True,
                    cuff_follow=True,
                    sort_bam=True
                    ),
                    bash.mkdir(link_directory),
                    bash.ln(
                        os.path.relpath(sample_log, link_directory),
                        os.path.join(link_directory, readset.sample.name + "_star.Log.final.out"),
                        sample_log
                    )
                ]
            )
            self.multiqc_inputs.append(sample_log)

            job.samples = [readset.sample]
            job.readsets = [readset]
            if not mapping or mapping == "2-pass":
                job.input_files.append(os.path.join(project_index_directory, "SAindex"))

            # If this readset is unique for this sample, further BAM merging is not necessary.
            # Thus, create a sample BAM symlink to the readset BAM.
            # remove older symlink before otherwise it raise an error if the link already exist (in case of redo)
            if len(readset.sample.readsets) == 1:
                readset_bam = os.path.join(alignment_pass_directory, "Aligned.sortedByCoord.out.bam")
                sample_bam = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name, readset.sample.name + ".sorted.bam")

                job = concat_jobs(
                    [
                        job,
                        bash.ln(
                            os.path.relpath(readset_bam, os.path.dirname(sample_bam)),
                            sample_bam,
                            input=readset_bam
                        )
                    ],
                    removable_files=[sample_bam]
                )

            self.multiqc_inputs.append(sample_log)
            job.name = "star_align.2." + readset.name if not mapping or mapping == "2-pass" else "star_align." + readset.name
            job.samples = [readset.sample]
            job.readsets = [readset]
            jobs.append(job)

        return jobs

    def picard_merge_sam_files(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        for sample in self.samples:
            # Skip samples with one readset only, since symlink has been created at align step
            if len(sample.readsets) > 1:
                alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
                inputs = [os.path.join(alignment_directory, readset.name, "Aligned.sortedByCoord.out.bam") for readset in sample.readsets]
                output = os.path.join(alignment_directory, sample.name + ".sorted.bam")

                jobs.append(
                    concat_jobs(
                        [
                            bash.rm(output),
                            picard.merge_sam_files(
                                inputs, 
                                output)
                        ],
                        name = "picard_merge_sam_files." + sample.name,
                        samples = [sample],
                        readsets = list(sample.readsets),
                        input_dependency = inputs
                    )
                )

        return jobs

    def picard_sort_sam(self):
        """
        The alignment file is reordered (QueryName) using [Picard](http://broadinstitute.github.io/picard/). The QueryName-sorted bam files will be used to determine raw read counts.
        """

        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name)

            job = picard.sort_sam(
                alignment_file_prefix + ".sorted.bam",
                alignment_file_prefix + ".QueryNameSorted.bam",
                "queryname"
            )
            job.name = "picard_sort_sam." + sample.name
            job.samples = [sample]
            jobs.append(job)
        return jobs

    def mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")
        for sample in self.samples:
            alignment_file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(link_directory),
                        picard.mark_duplicates(
                            [alignment_file_prefix + "bam"],
                            alignment_file_prefix + "mdup.bam",
                            alignment_file_prefix + "mdup.metrics",
                            create_index=False,
                            ini_section='mark_duplicates'
                        ),
                        sambamba.index(
                            alignment_file_prefix + "mdup.bam",
                            alignment_file_prefix + "mdup.bam.bai"
                        ),
                        bash.ln(
                            os.path.relpath(alignment_file_prefix + "mdup.metrics", link_directory),
                            os.path.join(link_directory, sample.name + ".sorted.mdup.metrics" ),
                            alignment_file_prefix + "mdup.metrics"
                            )
                    ],
                    name="mark_duplicates." + sample.name,
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
            )

            self.multiqc_inputs.append(alignment_file_prefix + "mdup.metrics")

        return jobs

    def picard_rna_metrics(self):
        """
        Computes a series of quality control metrics using both CollectRnaSeqMetrics and CollectAlignmentSummaryMetrics functions
        metrics are collected using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        reference_file = global_conf.global_get('picard_rna_metrics', 'genome_fasta', param_type='filepath')
        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")
        for sample in self.samples:
            alignment_file = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.mdup.bam")
            output_directory = os.path.join(self.output_dirs["metrics_directory"], sample.name)

            job_name = f"picard_rna_metrics.{sample.name}"
            job_project_tracking_metrics = []
            if self.project_tracking_json:
                job_project_tracking_metrics = concat_jobs(
                    [
                    picard.parse_aligned_reads_ratio_metrics_pt(os.path.join(output_directory, f"{sample.name}.alignment_summary_metrics")),
                    job2json_project_tracking.run(
                        input_file=os.path.join(output_directory, f"{sample.name}.alignment_summary_metrics"),
                        pipeline=self,
                        samples=sample.name,
                        readsets=",".join([readset.name for readset in sample.readsets]),
                        job_name=job_name,
                        metrics="aligned_reads_ratio=$aligned_reads_ratio"
                        )
                    ])

            job = concat_jobs(
                [
                    bash.mkdir(output_directory),
                    bash.mkdir(link_directory),
                    picard.collect_multiple_metrics(
                        alignment_file,
                        os.path.join(output_directory, sample.name),
                        reference_file,
                        library_type=sample.readsets[0].run_type
                    ),
                    picard.collect_rna_metrics(
                        alignment_file,
                        os.path.join(output_directory, sample.name+".picard_rna_metrics")
                    ),
                    bash.ln(
                        os.path.relpath(os.path.join(output_directory, sample.name + ".alignment_summary_metrics"), link_directory),
                        os.path.join(link_directory, sample.name + ".alignment_summary_metrics"),
                        os.path.join(output_directory, sample.name + ".alignment_summary_metrics")
                        ),
                    bash.ln(
                        os.path.relpath(os.path.join(output_directory, sample.name + ".quality_by_cycle_metrics"), link_directory),
                        os.path.join(link_directory, sample.name + ".quality_by_cycle_metrics"),
                        os.path.join(output_directory, sample.name + ".quality_by_cycle_metrics")
                        ),
                    bash.ln(
                        os.path.relpath(os.path.join(output_directory, sample.name + ".quality_distribution_metrics"), link_directory),
                        os.path.join(link_directory, sample.name + ".quality_distribution_metrics"),
                        os.path.join(output_directory, sample.name + ".quality_distribution_metrics")
                        ),
                    bash.ln(
                        os.path.relpath(os.path.join(output_directory, sample.name + ".picard_rna_metrics"), link_directory),
                        os.path.join(link_directory, sample.name + ".picard_rna_metrics"),
                        os.path.join(output_directory, sample.name + ".picard_rna_metrics")
                        ),
                    job_project_tracking_metrics
                ]
            )


            if sample.readsets[0].run_type == "PAIRED_END":
                job = concat_jobs(
                        [
                            job,
                            bash.ln(
                                os.path.relpath(os.path.join(output_directory, sample.name + ".base_distribution_by_cycle_metrics"), link_directory),
                                os.path.join(link_directory, sample.name + ".base_distribution_by_cycle_metrics"),
                                os.path.join(output_directory, sample.name + ".base_distribution_by_cycle_metrics")
                                ),
                            bash.ln(
                                os.path.relpath(os.path.join(output_directory, sample.name + ".insert_size_metrics"), link_directory),
                                os.path.join(link_directory, sample.name + ".insert_size_metrics"),
                                os.path.join(output_directory, sample.name + ".insert_size_metrics")
                                )
                        ]
                    )

            job.name = job_name
            job.removable_files = [output_directory]
            job.samples = [sample]
            job.readsets = list(sample.readsets)
            jobs.append(job)

            self.multiqc_inputs.append(output_directory)
        return jobs

    def rseqc(self):
        """
        Computes a series of quality control metrics using both CollectRnaSeqMetrics and CollectAlignmentSummaryMetrics functions.
        Metrics are collected using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            [input] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.recal.bam")],
                    [os.path.join(alignment_directory, sample.name + ".sorted.mdup.bam")]
                ]
            )
            output_directory = os.path.join(self.output_dirs["metrics_directory"], "rna", sample.name)

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            output_directory,
                            remove=True
                        ),
                        bash.chdir(
                            output_directory
                        ),
                        rseqc.tin(
                            os.path.relpath(input, output_directory),
                            output_directory
                        )
                    ],
                    name="rseqc.tin." + sample.name,
                    samples=[sample],
                    readsets=list(sample.readsets),
                    input_dependency=[input]
                )
            )

        return jobs


    def gatk_callable_loci(self):
        """
        Computes the callable region or the genome as a bed track.
        """

        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".")
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            [input] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.recal.bam")],
                    [os.path.join(alignment_directory, sample.name + ".sorted.mdup.bam")]
                ]
            )

            job = gatk4.callable_loci(
                input,
                alignment_file_prefix + "callable.bed",
                alignment_file_prefix + "callable.summary.txt"
            )
            job.name = "gatk_callable_loci." + sample.name
            job.samples = [sample]
            job.readsets = list(sample.readsets)
            jobs.append(job)

        return jobs

### THIS IS DEPRECATED. Replaced by rnaseqc2, should we remove it?
    def rnaseqc(self):
        """
        Computes a series of quality control metrics using [RNA-SeQC](https://www.broadinstitute.org/cancer/cga/rna-seqc).
        """

        jobs = []
        sample_file = os.path.join(self.output_dirs["alignment_directory"], "rnaseqc.samples.txt")
        sample_rows = [[sample.name, os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.mdup.bam"), "RNAseq"] for sample in self.samples]
        input_bams = [sample_row[1] for sample_row in sample_rows]
        output_directory = os.path.join(self.output_dirs["metrics_directory"], "rnaseqRep")
        # Use GTF with transcript_id only otherwise RNASeQC fails
        gtf_transcript_id = global_conf.global_get('rnaseqc', 'gtf_transcript_id', param_type='filepath')

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(output_directory),
                    Job(
                        input_bams,
                        [sample_file],
                        command="""\
echo "Sample\tBamFile\tNote
{sample_rows}" \\
  > {sample_file}""".format(
                            sample_rows="\n".join(["\t".join(sample_row) for sample_row in sample_rows]),
                            sample_file=sample_file
                        )
                    ),
                    metrics.rnaseqc(
                        sample_file,
                        output_directory,
                        self.run_type == "SINGLE_END",
                        gtf_file=gtf_transcript_id,
                        reference=global_conf.global_get('rnaseqc', 'genome_fasta', param_type='filepath'),
                        ribosomal_interval_file=global_conf.global_get('rnaseqc', 'ribosomal_fasta', param_type='filepath')
                    ),
                    bash.zip(
                        output_directory,
                        output_directory + ".zip",
                        recursive=True
                    )
                ],
                name="rnaseqc",
                removable_files=[output_directory],
                samples=self.samples,
                readsets=self.readsets
            )
        )

        trim_metrics_file = os.path.join(self.output_dirs["metrics_directory"], "trimSampleTable.tsv")
        metrics_file = os.path.join(self.output_dirs["metrics_directory"], "rnaseqRep", "metrics.tsv")
        report_metrics_file = os.path.join(self.output_dirs["report_directory"], "trimAlignmentTable.tsv")
        report_file = os.path.join(self.output_dirs["report_directory"], "RnaSeq.rnaseqc.md")
        jobs.append(
            Job(
                [metrics_file],
                [report_file, report_metrics_file],
                [
                    ['rnaseqc', 'module_python'],
                    ['rnaseqc', 'module_pandoc']],
                # Ugly awk to merge sample metrics with trim metrics if they exist; knitr may do this better
                command="""\
mkdir -p {report_dir} && \\
cp {output_directory}.zip {report_dir}/reportRNAseqQC.zip && \\
python -c 'import csv; csv_in = csv.DictReader(open("{metrics_file}"), delimiter="\t")
print "\t".join(["Sample", "Aligned Reads", "Alternative Alignments", "%", "rRNA Reads", "Coverage", "Exonic Rate", "Genes"])
print "\\n".join(["\t".join([
    line["Sample"],
    line["Mapped"],
    line["Alternative Aligments"],
    str(float(line["Alternative Aligments"]) / float(line["Mapped"]) * 100),
    line["rRNA"],
    line["Mean Per Base Cov."],
    line["Exonic Rate"],
    line["Genes Detected"]
]) for line in csv_in])' \\
  > {report_metrics_file}.tmp && \\
if [[ -f {trim_metrics_file} ]]
then
  awk -F"\t" 'FNR==NR{{raw_reads[$1]=$2; surviving_reads[$1]=$3; surviving_pct[$1]=$4; next}}{{OFS="\t"; if ($2=="Aligned Reads"){{surviving_pct[$1]="%"; aligned_pct="%"; rrna_pct="%"}} else {{aligned_pct=($2 / surviving_reads[$1] * 100); rrna_pct=($5 / surviving_reads[$1] * 100)}}; printf $1"\t"raw_reads[$1]"\t"surviving_reads[$1]"\t"surviving_pct[$1]"\t"$2"\t"aligned_pct"\t"$3"\t"$4"\t"$5"\t"rrna_pct; for (i = 6; i<= NF; i++) {{printf "\t"$i}}; print ""}}' \\
  {trim_metrics_file} \\
  {report_metrics_file}.tmp \\
  > {report_metrics_file}
else
  cp {report_metrics_file}.tmp {report_metrics_file}
fi && \\
rm {report_metrics_file}.tmp && \\
trim_alignment_table_md=`if [[ -f {trim_metrics_file} ]] ; then cut -f1-13 {report_metrics_file} | LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{print $1, sprintf("%\\47d", $2), sprintf("%\\47d", $3), sprintf("%.1f", $4), sprintf("%\\47d", $5), sprintf("%.1f", $6), sprintf("%\\47d", $7), sprintf("%.1f", $8), sprintf("%\\47d", $9), sprintf("%.1f", $10), sprintf("%.2f", $11), sprintf("%.2f", $12), sprintf("%\\47d", $13)}}}}' ; else cat {report_metrics_file} | LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{print $1, sprintf("%\\47d", $2), sprintf("%\\47d", $3), sprintf("%.1f", $4), sprintf("%\\47d", $5), sprintf("%.2f", $6), sprintf("%.2f", $7), $8}}}}' ; fi`
pandoc \\
  {report_template_dir}/{basename_report_file} \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable trim_alignment_table="$trim_alignment_table_md" \\
  --to markdown \\
  > {report_file}""".format(
                    output_directory=output_directory,
                    report_template_dir=self.report_template_dir,
                    trim_metrics_file=trim_metrics_file,
                    metrics_file=metrics_file,
                    basename_report_file=os.path.basename(report_file),
                    report_dir=self.output_dirs["report_directory"],
                    report_metrics_file=report_metrics_file,
                    report_file=report_file
                ),
                report_files=[report_file],
                name="rnaseqc_report",
                samples=self.samples,
                readsets=self.readsets
            )
        )
        return jobs

    def rnaseqc2(self):
        """
        Computes a series of quality control metrics using [RNA-SeQC](https://github.com/getzlab/rnaseqc).
        """
        jobs = []

        output_directory = os.path.join(self.output_dirs["metrics_directory"], "rna")
        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")
        scientific_name = global_conf.global_get('star', 'scientific_name')

        if os.path.isfile(os.path.join(global_conf.global_get('DEFAULT', 'annotations_prefix', param_type='prefixpath', required=False) + ".rnaseqc2.gtf")):
            input_gtf = os.path.join(global_conf.global_get('DEFAULT', 'annotations_prefix', param_type='prefixpath', required=False) + ".rnaseqc2.gtf")
        else:
            jobs.append(
                concat_jobs(
                    [
                    bash.mkdir(output_directory),
                    gtex_pipeline.collapse_gtf(
                        input_gtf=global_conf.global_get('rnaseqc2', 'gtf', param_type='filepath', required=True),
                        output_gtf=os.path.join(output_directory, scientific_name + "_collapsed.gtf")
                    )
                    ],
                    name="collapse_gtf"
                )
            )
            input_gtf = os.path.join(output_directory, scientific_name + "_collapsed.gtf")

        # Use GTF with transcript_id only otherwise RNASeQC fails
        for sample in self.samples:
            alignment_file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name)
            [input] = self.select_input_files(
                [
                    [os.path.join(alignment_file_prefix + ".sorted.mdup.bam")]
                ]
            )

            output_directory = os.path.join(self.output_dirs["metrics_directory"], "rna", sample.name, "rnaseqc2")
            job_name = f"rnaseqc2.{sample.name}"
            job_project_tracking_metrics = []
            if self.project_tracking_json:
                job_project_tracking_metrics = concat_jobs(
                    [
                    rnaseqc2.parse_expression_profiling_efficiency_metrics_pt(os.path.join(output_directory, sample.name + ".sorted.mdup.bam.metrics.tsv")),
                    job2json_project_tracking.run(
                        input_file=os.path.join(output_directory, sample.name + ".sorted.mdup.bam.metrics.tsv"),
                        pipeline=self,
                        samples=sample.name,
                        readsets=",".join([readset.name for readset in sample.readsets]),
                        job_name=job_name,
                        metrics="expression_profiling_efficiency=$expression_profiling_efficiency"
                        ),
                    rnaseqc2.parse_rrna_rate_metrics_pt(os.path.join(output_directory, sample.name + ".sorted.mdup.bam.metrics.tsv")),
                    job2json_project_tracking.run(
                        input_file=os.path.join(output_directory, sample.name + ".sorted.mdup.bam.metrics.tsv"),
                        pipeline=self,
                        samples=sample.name,
                        readsets=",".join([readset.name for readset in sample.readsets]),
                        job_name=job_name,
                        metrics="rrna_rate=$rrna_rate"
                        )
                    ])
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(output_directory),
                        bash.mkdir(link_directory),
                        rnaseqc2.run(
                            input,
                            input_gtf,
                            output_directory,
                        ),
                        bash.ln(
                            os.path.relpath(os.path.join(output_directory, sample.name + ".sorted.mdup.bam.metrics.tsv"), link_directory),
                            os.path.join(link_directory, sample.name + ".sorted.mdup.bam.metrics.tsv"),
                            os.path.join(output_directory, sample.name + ".sorted.mdup.bam.metrics.tsv")
                        ),
                        job_project_tracking_metrics
                    ],
                    name=job_name,
                    removable_files=[output_directory],
                    samples=self.samples,
                    readsets=self.readsets
                )
            )
            self.multiqc_inputs.append(os.path.join(output_directory, sample.name + ".sorted.mdup.bam.metrics.tsv"))
        return jobs

    def split_N_trim(self):
        """
        SplitNtrim. A [GATK](https://software.broadinstitute.org/gatk/) tool called SplitNCigarReads developed specially for RNAseq, which splits reads into exon segments (getting rid of Ns but maintaining grouping information) and hard-clip any sequences overhanging into the intronic regions.
        """

        jobs = []

        nb_jobs = global_conf.global_get('gatk_split_N_trim', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of haplotype jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for sample in self.samples:
            alignment_dir = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            split_dir = os.path.join(alignment_dir, "splitNtrim")
            split_file_prefix = os.path.join(split_dir, sample.name + ".")
            alignment_file_prefix = os.path.join(alignment_dir, sample.name + ".")

            if nb_jobs == 1:
                job = concat_jobs(
                    [
                        bash.mkdir(
                            split_dir,
                            remove=True
                        ),
                        gatk4.split_n_cigar_reads(
                            alignment_file_prefix + "sorted.mdup.bam",
                            split_file_prefix + "sorted.mdup.split.bam"
                        ),
                        bash.ln(
                            os.path.relpath(split_file_prefix + "sorted.mdup.split.bam", alignment_dir),
                            alignment_file_prefix + "sorted.mdup.split.bam",
                            input=split_file_prefix + "sorted.mdup.split.bam"
                        )
                    ],
                    name="gatk_split_N_trim." + sample.name,
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
                jobs.append(job)

            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):
                    job = concat_jobs(
                        [
                            # Create output directory since it is not done by default by GATK tools
                            bash.mkdir(
                                split_dir,
                                remove=True
                            ),
                            gatk4.split_n_cigar_reads(
                                alignment_file_prefix + "sorted.mdup.bam",
                                split_file_prefix + "sorted.mdup.split." + str(idx) + ".bam",
                                intervals=sequences
                            )
                        ],
                        name="gatk_split_N_trim." + sample.name + "." + str(idx),
                        samples=[sample],
                        readsets=list(sample.readsets)
                    )
                    jobs.append(job)

                job = concat_jobs(
                    [
                        bash.mkdir(
                            split_dir,
                            remove=True
                        ),
                        gatk4.split_n_cigar_reads(
                            alignment_file_prefix + "sorted.mdup.bam",
                            split_file_prefix + "sorted.mdup.split.others.bam",
                            exclude_intervals=unique_sequences_per_job_others
                        )
                    ],
                    name="gatk_split_N_trim." + sample.name + ".others",
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
                jobs.append(job)

        return jobs

    def sambamba_merge_splitNtrim_files(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [Sambamba](http://lomereiter.github.io/sambamba/docs/sambamba-merge.html).
        """

        jobs = []

        nb_jobs = global_conf.global_get('gatk_split_N_trim', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of haplotype jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            split_file_prefix = os.path.join(alignment_directory, "splitNtrim", sample.name + ".")
            output = os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.bam")

            if nb_jobs > 1:
                unique_sequences_per_job, _ = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                inputs = []
                for idx, _ in enumerate(unique_sequences_per_job):
                    inputs.append(split_file_prefix + "sorted.mdup.split." + str(idx) + ".bam")
                inputs.append(split_file_prefix + "sorted.mdup.split.others.bam")

                job = sambamba.merge(inputs, output)
                job.name = "sambamba_merge_splitNtrim_files." + sample.name
                job.samples = [sample]
                job.readsets = list(sample.readsets)
                jobs.append(job)

        return jobs

    def gatk_indel_realigner(self):
        """
        Insertion and deletion realignment is performed on regions where multiple base mismatches
        are preferred over indels by the aligner since it can appear to be less costly by the algorithm.
        Such regions will introduce false positive variant calls which may be filtered out by realigning
        those regions properly. Realignment is done using [GATK](https://www.broadinstitute.org/gatk/).
        The reference genome is divided by a number regions given by the `nb_jobs` parameter.
        """

        jobs = []

        nb_jobs = global_conf.global_get('gatk_indel_realigner', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of realign jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            input = os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.bam")

            quality_offsets = [readset.quality_offset for readset in sample.readsets]

            if nb_jobs == 1:
                realign_prefix = os.path.join(realign_directory, "all")
                realign_intervals = realign_prefix + ".intervals"
                output_bam = realign_prefix + ".bam"
                sample_output_bam = os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.bam")
                job = concat_jobs(
                    [
                        bash.mkdir(
                            realign_directory,
                            remove=True
                        ),
                        gatk4.realigner_target_creator(
                            input,
                            realign_intervals,
                            output_dir=self.output_dir,
                            fix_encoding=True if quality_offsets[0] == 64 else ""
                        ),
                        gatk4.indel_realigner(
                            input,
                            output=output_bam,
                            target_intervals=realign_intervals,
                            output_dir=self.output_dir,
                            fix_encoding=True if quality_offsets[0] == 64 else ""
                        ),
                        # Create sample realign symlink since no merging is required
                        bash.ln(
                            os.path.relpath(output_bam, alignment_directory),
                            sample_output_bam,
                            input=output_bam
                        )
                    ],
                    name="gatk_indel_realigner." + sample.name,
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
                jobs.append(job)

            else:
                # The first sequences are the longest to process.
                # Each of them must be processed in a separate job.
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                # Create one separate job for each of the first sequences
                for idx, sequence in enumerate(unique_sequences_per_job):
                    realign_prefix = os.path.join(realign_directory, str(idx))
                    realign_intervals = realign_prefix + ".intervals"
                    intervals = sequence
                    if str(idx) == 0:
                        intervals.append("unmapped")
                    output_bam = realign_prefix + ".bam"
                    job = concat_jobs(
                        [
                            # Create output directory since it is not done by default by GATK tools
                            bash.mkdir(
                                realign_directory,
                                remove=True
                            ),
                            gatk4.realigner_target_creator(
                                input,
                                realign_intervals,
                                intervals=intervals,
                                output_dir=self.output_dir,
                                fix_encoding=True if quality_offsets[0] == 64 else ""
                            ),
                            gatk4.indel_realigner(
                                input,
                                output=output_bam,
                                target_intervals=realign_intervals,
                                intervals=intervals,
                                output_dir=self.output_dir,
                                fix_encoding=True if quality_offsets[0] == 64 else ""
                            )
                        ],
                        name="gatk_indel_realigner." + sample.name + "." + str(idx),
                        samples=[sample],
                        readsets=list(sample.readsets)
                    )
                    jobs.append(job)

                # Create one last job to process the last remaining sequences and 'others' sequences
                realign_prefix = os.path.join(realign_directory, "others")
                realign_intervals = realign_prefix + ".intervals"
                output_bam = realign_prefix + ".bam"
                job = concat_jobs(
                    [
                        # Create output directory since it is not done by default by GATK tools
                        bash.mkdir(
                            realign_directory,
                            remove=True
                        ),
                        gatk4.realigner_target_creator(
                            input,
                            realign_intervals,
                            exclude_intervals=unique_sequences_per_job_others,
                            output_dir=self.output_dir,
                            fix_encoding=True if quality_offsets[0] == 64 else ""
                        ),
                        gatk4.indel_realigner(
                            input,
                            output=output_bam,
                            target_intervals=realign_intervals,
                            exclude_intervals=unique_sequences_per_job_others,
                            output_dir=self.output_dir,
                            fix_encoding=True if quality_offsets[0] == 64 else ""
                        )
                    ],
                    name="gatk_indel_realigner." + sample.name + ".others",
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
                jobs.append(job)

        return jobs

    def sambamba_merge_realigned(self):
        """
        BAM files of regions of realigned reads are merged per sample using [Sambamba](http://lomereiter.github.io/sambamba/index.html).
        """

        jobs = []

        nb_jobs = global_conf.global_get('gatk_indel_realigner', 'nb_jobs', param_type='posint')

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            merged_realigned_bam = os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.bam")

            # if nb_jobs == 1, symlink has been created in indel_realigner and merging is not necessary
            if nb_jobs > 1:
                unique_sequences_per_job, _ = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                inputs = []
                for idx, _ in enumerate(unique_sequences_per_job):
                    inputs.append(os.path.join(realign_directory, str(idx) + ".bam"))
                inputs.append(os.path.join(realign_directory, "others.bam"))

                job = sambamba.merge(inputs, merged_realigned_bam)
                job.name = "sambamba_merge_realigned." + sample.name
                job.removable_files = [merged_realigned_bam]
                job.samples = [sample]
                job.readsets = list(sample.readsets)
                jobs.append(job)

        return jobs

    def recalibration(self):
        """
        Recalibrate base quality scores of sequencing-by-synthesis reads in an aligned BAM file. After recalibration,
        the quality scores in the QUAL field in each read in the output BAM are more accurate in that
        the reported quality score is closer to its actual probability of mismatching the reference genome.
        Moreover, the recalibration tool attempts to correct for variation in quality with machine cycle
        and sequence context, and by doing so, provides not only more accurate quality scores but also
        more widely dispersed ones.
        """

        jobs = []

        created_interval_lists = []

        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            [input] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.bam")],
                    [os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.bam")]
                ]
            )
            print_reads_output = re.sub(r"\.bam$", ".recal.bam", input)
            base_recalibrator_output = re.sub(r"\.bam$", ".recalibration_report.grp", input)

            interval_list = None

            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            if coverage_bed:
                interval_list = re.sub(r"\.[^.]+$", ".interval_list", coverage_bed)

                if not interval_list in created_interval_lists:
                    job = tools.bed2interval_list(coverage_bed, interval_list)
                    job.name = "interval_list." + os.path.basename(coverage_bed)
                    job.samples = [sample]
                    job.readsets = list(sample.readsets)
                    jobs.append(job)
                    # created_interval_lists.append(interval_list)

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(link_directory),
                        gatk4.base_recalibrator(
                            input,
                            base_recalibrator_output,
                            intervals=interval_list
                        ),
                        bash.ln(
                            os.path.relpath(base_recalibrator_output, link_directory),
                            os.path.join(link_directory, sample.name + ".sorted.mdup.split.realigned.recalibration_report.grp"),
                            base_recalibrator_output
                            )
                    ],
                    name="gatk_base_recalibrator." + sample.name,
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        gatk4.print_reads(
                            input,
                            print_reads_output,
                            base_recalibrator_output
                        ),
                        deliverables.md5sum(
                            print_reads_output,
                            print_reads_output + ".md5",
                            self.output_dir
                            )
                    ],
                    name="gatk_print_reads." + sample.name,
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
            )

            self.multiqc_inputs.append(base_recalibrator_output)

        return jobs

    def gatk_haplotype_caller(self):
        """
        GATK haplotype caller for snps and small indels.
        """

        jobs = []

        nb_haplotype_jobs = global_conf.global_get('gatk_haplotype_caller', 'nb_jobs', param_type='posint')
        if nb_haplotype_jobs > 50:
            log.warning(
                "Number of haplotype jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            haplotype_directory = os.path.join(alignment_directory, "rawHaplotypeCaller")
            [input] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.recal.bam")],
                    [os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.bam")],
                    [os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.recal.bam")],
                    [os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.bam")]
                ]
            )

            if nb_haplotype_jobs == 1:
                job = concat_jobs(
                    [
                        # Create output directory since it is not done by default by GATK tools
                        bash.mkdir(haplotype_directory),
                        gatk4.haplotype_caller(
                            input,
                            os.path.join(haplotype_directory, sample.name + ".hc.vcf.gz"),
                            interval_list=None
                        )
                    ],
                    removable_files=[haplotype_directory],
                    name="gatk_haplotype_caller." + sample.name,
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
                jobs.append(job)

            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(
                    self.sequence_dictionary,
                    nb_haplotype_jobs - 1,
                    variant=True
                )

                job = []

                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):
                    job = concat_jobs(
                        [
                            # Create output directory since it is not done by default by GATK tools
                            bash.mkdir(haplotype_directory),
                            gatk4.haplotype_caller(
                                input,
                                os.path.join(haplotype_directory, sample.name + "." + str(idx) + ".hc.vcf.gz"),
                                interval_list=None,
                                intervals=sequences
                            )
                        ],
                        removable_files=[haplotype_directory],
                        name="gatk_haplotype_caller." + sample.name + "." + str(idx),
                        samples=[sample],
                        readsets=list(sample.readsets)
                    )
                    jobs.append(job)

                # Create one last job to process the last remaining sequences and 'others' sequences
                job = concat_jobs(
                    [
                        # Create output directory since it is not done by default by GATK tools
                        bash.mkdir(haplotype_directory),
                        gatk4.haplotype_caller(
                            input,
                            os.path.join(haplotype_directory, sample.name + ".others.hc.vcf.gz"),
                            interval_list=None,
                            exclude_intervals=unique_sequences_per_job_others
                        )
                    ],
                    removable_files=[haplotype_directory],
                    name="gatk_haplotype_caller." + sample.name + ".others",
                    samples = [sample]
                )
                jobs.append(job)

        return jobs

    def merge_hc_vcf(self):
        """
        Merges vcfs from haplotype caller to generate a sample level vcf
        """

        jobs = []
        nb_haplotype_jobs = global_conf.global_get('gatk_haplotype_caller', 'nb_jobs', param_type='posint')

        for sample in self.samples:
            haplotype_file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, "rawHaplotypeCaller", sample.name)
            output_haplotype_file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name)
            if nb_haplotype_jobs == 1:
                vcfs_to_merge = [haplotype_file_prefix + ".hc.vcf.gz"]
            else:
                unique_sequences_per_job, _ = split_by_size(self.sequence_dictionary, nb_haplotype_jobs - 1)

                vcfs_to_merge = [haplotype_file_prefix + "." + str(idx) + ".hc.vcf.gz" for idx in range(len(unique_sequences_per_job))]
                vcfs_to_merge.append(haplotype_file_prefix + ".others.hc.vcf.gz")

            job = concat_jobs(
                [
                    pipe_jobs(
                        [
                            bcftools.concat(
                                vcfs_to_merge,
                                None
                            ),
                            bash.grep(
                                None,
                                None,
                                "-v 'GL00'"
                            ),
                            bash.grep(
                                None,
                                None,
                                "-Ev 'chrUn|random'"
                            ),
                            bash.grep(
                                None,
                                None,
                                "-vE 'EBV|hs37d5'"
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_haplotype_file_prefix + ".hc.vcf.gz"
                            )
                        ]
                    )
                ],
                name="merge_vcfs." + sample.name,
                samples=[sample],
                readsets=list(sample.readsets)
            )
            jobs.append(job)

        return jobs

    def run_vcfanno(self):
        """
        [vcfanno](https://github.com/brentp/vcfanno) is used to annotate VCF files with preferred INFO fields from anu number of VCF or BED files.
        """

        jobs = []

        for sample in self.samples:
            input = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".hc.vcf.gz")
            output = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".hc.rnaedit.vcf.gz")

            job = pipe_jobs(
                [
                    vcfanno.run(
                        input,
                        None,
                    ),
                    htslib.bgzip_tabix(
                        None,
                        output
                    )
                ],
                name="run_vcfanno.rnaedit." + sample.name,
                samples=[sample],
                readsets=list(sample.readsets)
            )
            jobs.append(job)

        return jobs

    def variant_filtration(self):
        """
        GATK VariantFiltration.
        VariantFiltration is a GATK tool for hard-filtering variant calls based on certain criteria. Records are hard-filtered
        by changing the value in the FILTER field to something other than PASS.
        """

        jobs = []

        for sample in self.samples:
            input_file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name)

            job = concat_jobs(
                [
                    gatk4.variant_filtration(
                        input_file_prefix + ".hc.rnaedit.vcf.gz",
                        input_file_prefix + ".hc.flt.vcf.gz",
                        global_conf.global_get('gatk_variant_filtration', 'other_options')
                    )
                ],
                name="gatk_variant_filtration." + sample.name,
                samples=[sample],
                readsets=list(sample.readsets)
            )
            jobs.append(job)

        return jobs

    def decompose_and_normalize(self):
        """
        [vt](https://genome.sph.umich.edu/wiki/Vt#Normalization) is used to normalized and decompose VCF files. For more
        information about normalizing and decomposing visit
        [here](https://research-help.genomicsengland.co.uk/display/GERE/Variant+Normalisation). An indexed file is also
        generated from the output file using [htslib](http://www.htslib.org/download/).
        """

        jobs = []

        for sample in self.samples:
            input_file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name)

            [input] = self.select_input_files(
                [
                    [os.path.join(input_file_prefix + ".hc.rnaedit.vcf.gz")],
                    [os.path.join(input_file_prefix + ".hc.flt.vcf.gz")]
                ]
            )

            job = pipe_jobs(
                [
                    vt.decompose_and_normalize_mnps(
                        input,
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        input_file_prefix + ".hc.vt.vcf.gz"
                    )
                ],
                name="decompose_and_normalize." + sample.name,
                samples=[sample],
                readsets=list(sample.readsets)
            )
            jobs.append(job)

        return jobs

    def compute_snp_effects(self):
        """
        [SnpEff](https://pcingola.github.io/SnpEff/) is used to variant annotation and effect prediction on genes by
        using an interval forest approach. It annotates and predicts the effects of genetic variants
        (such as amino acid changes).
        """
        jobs = []

        for sample in self.samples:
            input_file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name)

            job = concat_jobs(
                [
                    snpeff.compute_effects(
                        input_file_prefix + ".hc.vt.vcf.gz",
                        input_file_prefix + ".hc.snpeff.vcf",
                        options=global_conf.global_get('compute_effects', 'options')
                    ),
                    htslib.bgzip_tabix(
                        input_file_prefix + ".hc.snpeff.vcf",
                        input_file_prefix + ".hc.snpeff.vcf.gz"
                    )
                ],
                name="compute_effects." + sample.name,
                samples=[sample],
                readsets=list(sample.readsets)
            )
            jobs.append(job)

        return jobs

    def gemini_annotations(self):
        """
        [gemini](https://github.com/arq5x/gemini) (GEnome MINIng) is used to integrative exploration of genetic
        variation and genome annotations. For more information
        visit the [Gemini documentation] (https://gemini.readthedocs.io/en/latest/).
        """

        jobs = []
        gemini_module = global_conf.global_get("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        for sample in self.samples:
            out_dir = os.path.join(self.output_dirs["alignment_directory"])
            input_file_prefix = os.path.join(out_dir, sample.name, sample.name)

            job = concat_jobs(
                [
                    gemini.gemini_annotations(
                        input_file_prefix + ".hc.snpeff.vcf.gz",
                        input_file_prefix + ".gemini." + gemini_version + ".db"
                    )
                ],
                name="gemini_annotations." + sample.name,
                samples = [sample]
            )
            jobs.append(job)

        return jobs

    def filter_gatk(self):
        """
        Applies custom script to inject FORMAT information - tumor/normal DP and VAP into the INFO field
        the filter on those generated fields.
        """
        jobs = []

        for sample in self.samples:
            input = os.path.join(
                self.output_dirs["alignment_directory"],
                sample.name,
                sample.name + ".hc.vt.vcf.gz"
            )
            output_annot = os.path.join(
                self.output_dirs["alignment_directory"],
                sample.name,
                sample.name + ".hc.vt.annot.vcf.gz"
            )
            output_filter = os.path.join(
                self.output_dirs["alignment_directory"],
                sample.name,
                sample.name + ".hc.vt.annot.flt.vcf.gz"
            )

            jobs.append(
                concat_jobs(
                    [
                        tools.format2pcgr(
                            input,
                            output_annot,
                            global_conf.global_get('filter_gatk', 'call_filter'),
                            "somatic",
                            sample.name,
                            ini_section='filter_gatk'
                        ),
                        pipe_jobs(
                            [
                                bcftools.view(
                                    output_annot,
                                    None,
                                    filter_options=global_conf.global_get('filter_gatk', 'somatic_filter_options'),
                                ),
                                bcftools.sort(
                                    None,
                                    output_filter,
                                    sort_options="-Oz"
                                )
                            ]
                        ),
                        htslib.tabix(
                            output_filter,
                            options="-pvcf"
                        ),
                        deliverables.md5sum(
                            output_filter,
                            output_filter + ".md5",
                            self.output_dir
                            )
                    ],
                    name="filter_gatk." + sample.name,
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
            )

        return jobs

    def report_cpsr(self):
        """
        Creates a cpsr gremline report (https://sigven.github.io/cpsr/)
        input: filtered ensemble gremline vcf
        output: html report and addtionalflat files
        """
        jobs = []

        for sample in self.samples:
            input = os.path.join(
                self.output_dirs["alignment_directory"],
                sample.name,
                sample.name + ".hc.vt.annot.flt.vcf.gz"
            )

            cpsr_directory = os.path.join(
                self.output_dirs["alignment_directory"],
                sample.name,
                "cpsr"
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(cpsr_directory),
                        cpsr.report(
                            input,
                            cpsr_directory,
                            sample.name
                        )
                    ],
                    name="report_cpsr." + sample.name,
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
            )

        return jobs

    def report_pcgr(self):
        """
        Creates a PCGR somatic + germline report (https://sigven.github.io/cpsr/)
        input: filtered ensemble gremline vcf
        output: html report and addtionalflat files
        """
        jobs = []

        assembly = global_conf.global_get('report_pcgr', 'assembly')

        for sample in self.samples:
            cpsr_directory = os.path.join(
                self.output_dirs["alignment_directory"],
                sample.name,
                "cpsr"
            )
            input_cpsr = os.path.join(
                cpsr_directory,
                sample.name + ".cpsr." + assembly + ".json.gz"
            )
            input = os.path.join(
                self.output_dirs["alignment_directory"],
                sample.name,
                sample.name + ".hc.vt.annot.flt.vcf.gz"
            )
            pcgr_directory = os.path.join(
                self.output_dirs["alignment_directory"],
                sample.name,
                "pcgr"
            )
            # PCGR (before v2) does not accept sample IDs longer than 35 characters and uses the sample ID to name output files.
            # For samples that have longer sample IDs the output files will have non-matching names, so create symlinks with full-length names.
            if global_conf.global_get('report_pcgr', 'module_pcgr').split("/")[2] < "2" and sample.name != sample.name[:35]:
                output = [
                        os.path.join(pcgr_directory, sample.name[:35] + ".pcgr_acmg." + assembly + ".flexdb.html"),
                        os.path.join(pcgr_directory, sample.name[:35] + ".pcgr_acmg." + assembly + ".maf"),
                        os.path.join(pcgr_directory, sample.name[:35] + ".pcgr_acmg." + assembly + ".snvs_indels.tiers.tsv")
                    ]
                final_command = concat_jobs(
                        [
                            bash.ln(
                                os.path.relpath(output[0], pcgr_directory),
                                os.path.join(pcgr_directory, sample.name + ".pcgr_acmg." + assembly + ".flexdb.html"),
                                input=output[0]
                                ),
                            bash.ln(
                                os.path.relpath(output[1], pcgr_directory),
                                os.path.join(pcgr_directory, sample.name + ".pcgr_acmg." + assembly + ".maf"),
                                input=output[1]
                                ),
                            bash.ln(
                                os.path.relpath(output[2], pcgr_directory),
                                os.path.join(pcgr_directory, sample.name + ".pcgr_acmg." + assembly + ".snvs_indels.tiers.tsv"),
                                input=output[2]
                                )
                            ]
                        )
            # output file name patterns have changed in pcgr versions >2.0.0
            # cpsr input changed in pcgr 2.1.0
            elif global_conf.global_get('report_pcgr', 'module_pcgr').split("/")[2] > "2":
                input_cpsr = os.path.join(
                    cpsr_directory,
                    sample.name + ".cpsr." + assembly
                )
                output = [
                        os.path.join(pcgr_directory, sample.name + ".pcgr." + assembly + ".html"),
                        os.path.join(pcgr_directory, sample.name + ".pcgr." + assembly + ".maf"),
                        os.path.join(pcgr_directory, sample.name + ".pcgr." + assembly + ".snvs_indels_ann.tsv.gz")
                    ]
                final_command = bash.ls(output[0])
                input_dependencies = [input, input_cpsr + ".classification.tsv.gz", input_cpsr + ".conf.yaml"]
            
            else:
                output = [
                        os.path.join(pcgr_directory, sample.name + ".pcgr_acmg." + assembly + ".flexdb.html"),
                        os.path.join(pcgr_directory, sample.name + ".pcgr_acmg." + assembly + ".maf"),
                        os.path.join(pcgr_directory, sample.name + ".pcgr_acmg." + assembly + ".snvs_indels.tiers.tsv")
                    ]
                final_command = bash.ls(output[0])
                input_dependencies = [input, input_cpsr]

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            pcgr_directory,
                        ),
                        pcgr.report(
                            input,
                            input_cpsr,
                            pcgr_directory,
                            sample.name,
                            input_cna=None
                        ),
                        final_command
                    ],
                    name="report_pcgr." + sample.name,
                    samples=[sample],
                    readsets=list(sample.readsets),
                    input_dependency=input_dependencies,
                    output_dependency=output
                )
            )

        return jobs

    def run_star_fusion(self):
        """
        [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) is a component of the Trinity Cancer Transcriptome Analysis Toolkit (CTAT). Based on the STAR
        aligner it identifies candidate fusion transcripts supported by Illumina reads.
        """

        jobs = []

        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        left_fastqs = collections.defaultdict(list)
        right_fastqs = collections.defaultdict(list)

        for readset in self.readsets:
            trim_file_prefix = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + "-trimmed-")

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([
                        re.sub(r"\.bam$", ".pair1.fastq.gz", readset.bam),
                        re.sub(r"\.bam$", ".pair2.fastq.gz", readset.bam)
                        ])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub(r"\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None

            else:
                raise Exception("Error: run type \"" + readset.run_type + "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            left_fastqs[readset.sample.name].append(fastq1)
            right_fastqs[readset.sample.name].append(fastq2)

        for sample in self.samples:
            output_dir = os.path.join(self.output_dirs["fusion_directory"], sample.name, "star_fusion")

            star_fusion_job = concat_jobs(
                [
                    bash.mkdir(output_dir),
                    bash.mkdir(link_directory),
                    star_fusion.run(
                        left_fastqs[sample.name],
                        right_fastqs[sample.name],
                        output_dir
                    ),
                    bash.ln(
                        os.path.relpath(os.path.join(output_dir, "Log.final.out"), link_directory),
                        os.path.join(link_directory, sample.name + "_star_fusion.Log.final.out"),
                        os.path.join(output_dir, "Log.final.out")
                        )
                ]
                )

            if global_conf.global_get('run_star_fusion', 'force') == "True":
                job = concat_jobs(
                        [
                            bash.rm(
                                os.path.join(output_dir, "_starF_checkpoints")
                            ),
                            star_fusion_job
                        ],
                        input_dependency=left_fastqs[sample.name] + right_fastqs[sample.name]
                    )
            else:
                job = star_fusion_job

            job.name = f"run_star_fusion.{sample.name}"
            job.samples = [sample]
            job.readsets = list(sample.readsets)
            jobs.append(job)

            self.multiqc_inputs.append(os.path.join(output_dir, "Log.final.out"))
        return jobs

    def run_arriba(self):
        """
        [arriba](https://github.com/suhrig/arriba) is used for the detection of gene fusions from RNA-Seq data.
        arriba is based on the [STAR](https://github.com/alexdobin/STAR) aligner. Apart from gene fusions,
        Arriba can detect other structural rearrangements with potential clinical relevance, including viral integration
        sites, internal tandem duplications, whole exon duplications and intragenic inversions etc...
        """

        jobs = []

        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        left_fastqs = collections.defaultdict(list)
        right_fastqs = collections.defaultdict(list)

        for readset in self.readsets:
            trim_file_prefix = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + "-trimmed-")

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append(
                        [
                            re.sub(r"\.bam$", ".pair1.fastq.gz", readset.bam),
                            re.sub(r"\.bam$", ".pair2.fastq.gz", readset.bam)
                        ]
                    )
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub(r"\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None

            else:
                raise Exception("Error: run type \"" + readset.run_type + "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            left_fastqs[readset.sample.name].append(fastq1)
            right_fastqs[readset.sample.name].append(fastq2)

        for sample in self.samples:
            output_dir = os.path.join(self.output_dirs["fusion_directory"], sample.name, "arriba")

            job = concat_jobs(
                [
                    bash.rm(output_dir),
                    bash.mkdir(output_dir),
                    bash.mkdir(link_directory),
                    bash.chdir(output_dir),
                    arriba.run(
                        [os.path.relpath(fastq1, output_dir) for fastq1 in left_fastqs[sample.name]],
                        [os.path.relpath(fastq2, output_dir) for fastq2 in right_fastqs[sample.name] if fastq2],
                        output_dir
                    ),
                    bash.chdir(self.output_dir),
                    bash.ln(
                        os.path.relpath(os.path.join(output_dir, "Log.final.out"), link_directory),
                        os.path.join(link_directory, sample.name + "_arriba.Log.final.out"),
                        os.path.join(output_dir, "Log.final.out")
                        )
                ],
                input_dependency=left_fastqs[sample.name] + right_fastqs[sample.name],
                output_dependency=[
                    os.path.join(output_dir, "Log.final.out"), os.path.join(output_dir, "fusions.tsv"),
                    os.path.join(link_directory, sample.name + "_arriba.Log.final.out"),
                    link_directory
                    ],
                name="run_arriba." + sample.name,
                samples=[sample],
                readsets=list(sample.readsets)
            )
            jobs.append(job)

            self.multiqc_inputs.append(os.path.join(output_dir, "Log.final.out"))
        return jobs

    def run_annofuse(self):
        """
        [annofuse](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03922-7)
        is a R package and it is used to annotate, prioritize, and interactively explore putative oncogenic
        RNA fusions.
        """
        jobs = []

        for sample in self.samples:
            arriba_output = os.path.join(self.output_dirs["fusion_directory"], sample.name, "arriba", "fusions.tsv")
            star_fusion_output = os.path.join(self.output_dirs["fusion_directory"], sample.name, "star_fusion", "star-fusion.fusion_predictions.abridged.coding_effect.tsv")
            output = os.path.join(self.output_dirs["fusion_directory"], sample.name, "annoFuse")

            job = concat_jobs(
                [
                    bash.mkdir(
                        output,
                        remove=True
                    ),
                    annoFuse.run(
                        arriba_output,
                        star_fusion_output,
                        os.path.join(output, sample.name)
                    )
                ],
                name="run_annoFuse." + sample.name,
                samples=[sample],
                readsets=list(sample.readsets)
            )

            job.samples = [sample]
            job.readsets = list(sample.readsets)
            jobs.append(job)

        return jobs

    def estimate_ribosomal_rna(self):
        """
        Use bwa mem to align reads on the rRNA reference fasta and count the number of read mapped
        The filtered reads are aligned to a reference fasta file of ribosomal sequence. The alignment is done per sequencing readset.
        The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
        BWA output BAM files are then sorted by coordinate using [Picard](http://broadinstitute.github.io/picard/).

        This step takes as input files: readset Bam files.
        """

        jobs = []

        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        for readset in self.readsets:
            readset_bam = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name, readset.name , "Aligned.sortedByCoord.out.bam")
            output_folder = os.path.join(self.output_dirs["metrics_directory"], readset.sample.name, readset.name)
            readset_metrics_bam = os.path.join(output_folder, readset.name +"rRNA.bam")

            job = concat_jobs(
                [
                    bash.mkdir(os.path.dirname(readset_bam)),
                    bash.mkdir(link_directory),
                    bash.mkdir(output_folder),
                    pipe_jobs(
                        [
                            bvatools.bam2fq(
                                readset_bam
                            ),
                            bwa.mem(
                                "/dev/stdin",
                                None,
                                read_group="'@RG" + \
                                    "\tID:" + readset.name + \
                                    "\tSM:" + readset.sample.name + \
                                    ("\tLB:" + readset.library if readset.library else "") + \
                                    ("\tPU:run" + readset.run + "_" + readset.lane if readset.run and readset.lane else "") + \
                                    ("\tCN:" + global_conf.global_get('bwa_mem_rRNA', 'sequencing_center') if global_conf.global_get('bwa_mem_rRNA', 'sequencing_center', required=False) else "") + \
                                    "\tPL:Illumina" + \
                                    "'",
                                ref=global_conf.global_get('bwa_mem_rRNA', 'ribosomal_fasta'),
                                ini_section='bwa_mem_rRNA'
                            ),
                            picard.sort_sam(
                                "/dev/stdin",
                                readset_metrics_bam,
                                "coordinate",
                                ini_section='picard_sort_sam_rrna'
                            )
                        ]
                    ),
                    tools.py_rrnaBAMcount(
                        bam=readset_metrics_bam,
                        gtf=global_conf.global_get('bwa_mem_rRNA', 'gtf'),
                        output=os.path.join(output_folder, readset.name+"rRNA.stats.tsv"),
                        typ="transcript"
                    ),
                    bash.ln(
                        target_file=os.path.relpath(os.path.join(output_folder, readset.name+"rRNA.stats.tsv"), link_directory),
                        link=os.path.join(link_directory,  readset.name+"rRNA.stats.tsv"),
                        input=os.path.join(output_folder, readset.name+"rRNA.stats.tsv")
                    )
                ],
                name=f"bwa_mem_rRNA.{readset.name}",
                samples=[readset.sample]
            )

            job.removable_files=[readset_metrics_bam]
            job.samples = [readset.sample]
            job.readsets = [readset]
            jobs.append(job)
            self.multiqc_inputs.append(os.path.join(output_folder, readset.name+"rRNA.stats.tsv"))
        return jobs

    def wiggle(self):
        """
        Generate wiggle tracks suitable for multiple browsers.
        """

        jobs = []

        for sample in self.samples:
            bam_file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.mdup.")
            input_bam = bam_file_prefix + "bam"
            tracks_dir = os.path.join(self.output_dirs["tracks_directory"])
            big_wig = os.path.join(self.output_dirs['tracks_directory'], "bigWig")

            if global_conf.global_get('wiggle', 'separate_strand') == 'NO':

                output_file = os.path.join(self.output_dirs["tracks_directory"], "bigWig", sample.name + ".bw")

                job= concat_jobs(
                    [
                        bash.mkdir(tracks_dir),
                        bash.mkdir(big_wig),
                        deeptools.bamcoverage(
                            input_bam,
                            output_file,

                        )
                    ],
                    name=f"wiggle.{sample.name}",
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
                jobs.append(job)

            else:

                output_file = os.path.join(self.output_dirs["tracks_directory"], "bigWig", sample.name + ".bw")

                job = concat_jobs(
                    [
                        bash.mkdir(tracks_dir),
                        bash.mkdir(big_wig),
                        deeptools.bamcoverage(
                            input_bam,
                            output_file,

                        )
                    ],
                    name="wiggle."+ sample.name,
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
                jobs.append(job)

                strand = "forward"
                output_file = os.path.join(self.output_dirs["tracks_directory"], "bigWig", f"{sample.name}_{strand}.bw")

                job = concat_jobs(
                    [
                        bash.mkdir(tracks_dir),
                        bash.mkdir(big_wig),
                        deeptools.bamcoverage(
                            input_bam,
                            output_file,
                            strand
                        )
                    ],
                    name="wiggle."+sample.name+"_"+strand,
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
                jobs.append(job)

                strand = "reverse"
                output_file = os.path.join(self.output_dirs["tracks_directory"], "bigWig", f"{sample.name}_{strand}.bw")

                job = concat_jobs(
                    [
                        bash.mkdir(tracks_dir),
                        bash.mkdir(big_wig),
                        deeptools.bamcoverage(
                            input_bam,
                            output_file,
                            strand
                        )
                    ],
                    name="wiggle."+sample.name+"_"+strand,
                    samples=[sample],
                    readsets=list(sample.readsets)
                )
                jobs.append(job)

        return jobs

    def raw_counts(self):
        """
        Count reads in features using [htseq-count](https://htseq.readthedocs.io/en/master/htseqcount.html).
        """

        jobs = []
        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        for sample in self.samples:
            alignment_file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name)
            input_bam = alignment_file_prefix + ".QueryNameSorted.bam"

            # Count reads
            output_count = os.path.join(self.output_dirs["raw_counts_directory"], sample.name + ".readcounts.tsv")
            stranded = global_conf.global_get('htseq_count', 'stranded')
            job = concat_jobs(
                [
                    bash.mkdir(f"{self.output_dirs['raw_counts_directory']}"),
                    bash.mkdir(link_directory),## ensure existance of link directory
                    htseq.htseq_count(
                        input_bam,
                        global_conf.global_get('htseq_count', 'gtf', param_type='filepath'),
                        output_count,
                        global_conf.global_get('htseq_count', 'options'),
                        stranded
                        ),
                    bash.ln(
                        target_file=os.path.relpath(output_count, link_directory),
                        link=os.path.join(link_directory, sample.name + ".tsv"),
                        input=output_count
                        )
                ],
                name="htseq_count." + sample.name,
                samples=[sample],
                readsets=list(sample.readsets)
            )
            jobs.append(job)

            self.multiqc_inputs.append(output_count)
        return jobs

    def raw_counts_metrics(self):
        """
        Create rawcount matrix, zip the wiggle tracks and create the saturation plots based on standardized read counts.
        """

        jobs = []

        # Create raw count matrix
        output_directory = self.output_dirs["DGE_directory"]
        read_count_files = [os.path.join(self.output_dirs["raw_counts_directory"], sample.name + ".readcounts.tsv") for sample in self.samples]
        output_matrix = os.path.join(output_directory, "rawCountMatrix.tsv")

        job = Job(
            read_count_files,
            [output_matrix],
            [['raw_counts_metrics', 'module_mugqic_tools']],
            name="metrics.matrix"
        )
        job.command = """\
mkdir -p {output_directory} && \\
gtf2tmpMatrix.awk \\
  {reference_gtf} \\
  {output_directory}/tmpMatrix.txt && \\
HEAD='Gene\\tSymbol' && \\
for read_count_file in \\
  {read_count_files}
do
  sort -k1,1 $read_count_file > {output_directory}/tmpSort.txt && \\
  join -1 1 -2 1 <(sort -k1,1 {output_directory}/tmpMatrix.txt) {output_directory}/tmpSort.txt > {output_directory}/tmpMatrix.2.txt && \\
  mv {output_directory}/tmpMatrix.2.txt {output_directory}/tmpMatrix.txt && \\
  na=$(basename $read_count_file | rev | cut -d. -f3- | rev) && \\
  HEAD="$HEAD\\t$na"
done && \\
echo -e $HEAD | cat - {output_directory}/tmpMatrix.txt | tr ' ' '\\t' > {output_matrix} && \\
rm {output_directory}/tmpSort.txt {output_directory}/tmpMatrix.txt""".format(
            reference_gtf=global_conf.global_get('raw_counts_metrics', 'gtf', param_type='filepath'),
            output_directory=output_directory,
            read_count_files=" \\\n  ".join(read_count_files),
            output_matrix=output_matrix
        )
        job.samples = self.samples
        job.readsets = self.readsets
        jobs.append(job)

        # Create Wiggle tracks archive

        wiggle_directory = os.path.join(self.output_dirs["tracks_directory"], "bigWig")
        wiggle_archive = os.path.join(self.output_dir, "tracks.zip")

        if global_conf.global_get('wiggle', 'separate_strand') == 'NO':
            wiggle_files = [os.path.join(wiggle_directory, sample.name + ".bw") for sample in self.samples]
        else:
            wiggle_files = [os.path.join(wiggle_directory, sample.name + ".bw") for sample in self.samples]
            strand="forward"
            wiggle_files = wiggle_files + [os.path.join(wiggle_directory, sample.name+"_"+strand +".bw") for sample in self.samples]
            strand="reverse"
            wiggle_files = wiggle_files + [os.path.join(wiggle_directory, sample.name+"_"+strand +".bw") for sample in self.samples]


        jobs.append(
            Job(
                wiggle_files,
                [wiggle_archive],
                name="metrics.wigzip",
                command=f"zip -r {wiggle_archive} {wiggle_directory}",
                samples=self.samples,
                readsets=self.readsets
                )
            )

        # RPKM and Saturation
        count_file = os.path.join(self.output_dirs["DGE_directory"], "rawCountMatrix.tsv")
        gene_size_file = global_conf.global_get('rpkm_saturation', 'gene_size', param_type='filepath')
        rpkm_directory = self.output_dirs["raw_counts_directory"]
        saturation_directory = os.path.join(self.output_dirs["metrics_directory"], "saturation")

        job = concat_jobs(
            [
                bash.mkdir(saturation_directory),
                bash.mkdir(rpkm_directory),
                metrics.rpkm_saturation(
                    count_file,
                    gene_size_file,
                    rpkm_directory,
                    saturation_directory
                )
            ],
            name="rpkm_saturation",
            samples=self.samples,
            readsets=self.readsets
        )
        jobs.append(job)

        return jobs


    def stringtie(self):
        """
        Assemble transcriptome using [stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml).
        """
        jobs = []

        gtf = global_conf.global_get('stringtie','gtf', param_type='filepath')

        for sample in self.samples:
            input_bam = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.mdup.bam")
            output_directory = os.path.join(self.output_dirs["stringtie_directory"], sample.name)

            job = stringtie.stringtie(input_bam, output_directory, gtf)
            job.name = "stringtie." + sample.name
            job.samples = [sample]
            job.readsets = list(sample.readsets)
            jobs.append(job)

        return jobs

    def stringtie_merge(self):
        """
        Merge assemblies into a master teranscriptome reference using [stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml).
        """

        jobs = []

        output_directory = os.path.join(self.output_dirs["stringtie_directory"], "AllSamples")
        sample_file = os.path.join(self.output_dirs["stringtie_directory"], "stringtie-merge.samples.txt")
        input_gtfs = [os.path.join(self.output_dirs["stringtie_directory"], sample.name, "transcripts.gtf") for sample in self.samples]
        gtf = global_conf.global_get('stringtie','gtf', param_type='filepath')

        if os.path.exists(os.path.join(self.output_dir, self.output_dirs["stringtie_directory"], "stringtieAbundDone")) and not self.force_jobs:
            log.info("Stringtie Abund done already... Skipping stringtie_merge step...")

        else:

            jobs = [
                concat_jobs(
                    [
                        bash.mkdir(output_directory),
                        Job(
                            input_gtfs,
                            [sample_file],
                            command="""\
`cat > {sample_file} << END
{sample_rows}
END

`""".format(
                                sample_rows="\n".join(input_gtfs),
                                sample_file=sample_file
                            )
                        ),
                        stringtie.stringtie_merge(
                            sample_file,
                            output_directory,
                            gtf
                        )
                    ],
                    name="stringtie_merge",
                    samples=self.samples,
                    readsets=self.readsets
                )
            ]

        return jobs

    def stringtie_abund(self):
        """
        Assemble transcriptome and compute RNA-seq expression using [stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml).
        """
        jobs = []

        gtf = os.path.join(self.output_dirs["stringtie_directory"], "AllSamples", "merged.gtf")

        donejob_input_dep = []
        for sample in self.samples:
            input_bam = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.mdup.bam")
            output_directory = os.path.join(self.output_dirs["stringtie_directory"], sample.name)

            job = stringtie.stringtie(
                input_bam,
                output_directory,
                gtf,
                abund=True
            )
            job.name = "stringtie_abund." + sample.name
            job.samples = [sample]
            job.readsets = list(sample.readsets)

            donejob_input_dep.extend(job.output_files)

            jobs.append(job)

        done_file = os.path.join(self.output_dirs["stringtie_directory"], "stringtieAbundDone")
        jobs.append(
            concat_jobs(
                [
                    bash.touch(done_file)
                ],
                name="stringtie_abund.donefile",
                samples=self.samples,
                readsets=self.readsets,
                input_dependency=donejob_input_dep,
                output_dependency=[done_file]
            )
        )
        return jobs

    def ballgown(self):
        """
        [Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html) is used to calculate differential transcript and gene expression levels and test them for significant differences.
        """

        jobs = []

        # Perform ballgown on each design contrast
        # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        if self.contrasts:
            design_file = os.path.relpath(self.design_file.name, self.output_dir)
        output_directory = self.output_dirs["ballgown_directory"]
        input_abund = [os.path.join(self.output_dirs["stringtie_directory"], sample.name, "abundance.tab") for sample in self.samples]

        ballgown_job = ballgown.ballgown(
            input_abund,
            design_file,
            output_directory,
            [os.path.join(output_directory, contrast.name, "gene_exp.diff") for contrast in self.contrasts] +
            [os.path.join(output_directory, contrast.name, "transcript_exp.diff") for contrast in self.contrasts]
        )
        ballgown_job.name = "ballgown"
        ballgown_job.samples = self.samples
        ballgown_job.readsets = self.readsets
        jobs.append(ballgown_job)

        return jobs

    def differential_expression(self):
        """
        Performs differential gene expression analysis using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html).
        Merge the results of the analysis in a single csv file.
        """

        # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        if self.contrasts:
            design_file = os.path.relpath(self.design_file.name, self.output_dir)

        output_directory = self.output_dirs["DGE_directory"]
        count_matrix = os.path.join(output_directory, "rawCountMatrix.tsv")
        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        edger_job = differential_expression.edger(
            design_file,
            count_matrix,
            output_directory,
            [os.path.join(output_directory, contrast.name, "edger_results.csv") for contrast in self.contrasts]
        )

        deseq_job = differential_expression.deseq2(
            design_file,
            count_matrix,
            output_directory,
            [os.path.join(output_directory, contrast.name, "dge_results.csv") for contrast in self.contrasts]
        )

        if self.batch_file:
            # If provided a batch file, compute DGE with batch effect correction
            batch_file = os.path.relpath(self.batch_file.name, self.output_dir)

            # edger_job_batch_corrected = differential_expression.edger(
            #     design_file,
            #     count_matrix,
            #     batch_file,
            #     [os.path.join(f"{output_directory}_batch_corrected", contrast.name, "edger_results.csv") for contrast in self.contrasts]
            #     f"{output_directory}_batch_corrected"
            # )

            deseq_job_batch_corrected = differential_expression.deseq2(
                design_file,
                count_matrix,
                f"{output_directory}_batch_corrected",
                [os.path.join(f"{output_directory}_batch_corrected", contrast.name, "deseq2_results.csv") for contrast in self.contrasts],
                batch_file
            )
        diff_exp_report_job = concat_jobs(
            [
                bash.mkdir(link_directory),
                tools.design_report(
                    design_file,
                    os.path.join(link_directory, "design_mqc.txt")
                )
            ]
        )
        self.multiqc_inputs.append(os.path.join(link_directory, "design_mqc.txt"))

        for contrast in self.contrasts:
            diff_exp_report_job = concat_jobs(
                [
                    diff_exp_report_job,
                    tools.diff_exp_report(
                        os.path.join(output_directory, contrast.name, "dge_results.csv"),
                        os.path.join(link_directory, contrast.name + ".dge_results_mqc.tsv"),
                        f"Differential Expression {contrast.name}",
                        contrast.name
                    )
                ]
            )
            self.multiqc_inputs.append(os.path.join(link_directory, contrast.name + ".dge_results_mqc.tsv"))

            

        return [
            concat_jobs(
                [
                    bash.mkdir(output_directory),
                    bash.mkdir(f"{output_directory}_batch_corrected") if self.batch_file else None,
                    edger_job,
                    deseq_job,
                    deseq_job_batch_corrected if self.batch_file else None,                    
                    diff_exp_report_job
                ],
                name="differential_expression",
                samples=self.samples,
                readsets=self.readsets
            )
        ]

### NO LONGER USED. Should we remove it?
    def ihec_metrics(self):
        """
        Generate IHEC's standard metrics.
        """

        genome = global_conf.global_get('ihec_metrics', 'assembly')
        return [
            metrics.ihec_metrics_rnaseq(
                [
                    os.path.join(self.output_dirs['metrics_directory'], "rnaseqRep", "metrics.tsv"),
                    os.path.join(self.output_dirs["report_directory"], "trimAlignmentTable.tsv")
                ] + [
                    os.path.join(self.output_dirs["metrics_directory"], readset.sample.name, readset.name, readset.name + "rRNA.stats.tsv")
                    for readset in self.readsets
                ],
                [os.path.join(self.output_dirs['report_directory'], "IHEC_metrics_rnaseq_All.txt")],
                genome
            )
        ]


    def multiqc(self):
        """
        Aggregate results from bioinformatics analyses across many samples into a single report.
        MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool,
        perfect for summarising the output from numerous bioinformatics tools [MultiQC](https://multiqc.info/).
        """
        jobs = []

        input_links = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")
        output = os.path.join(self.output_dirs['report_directory'], f"RnaSeq.{self.protocol}.multiqc")

        job = multiqc.run(
                [input_links],
                output
                )
        job.name = "multiqc"
        job.input_files = self.multiqc_inputs
        jobs.append(job)

        by_sample = global_conf.global_get("multiqc", "by_sample", required=False) # option to also generate individual multiqc reports for each sample

        if by_sample == "true":
            for sample in self.samples:
                input = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs", sample.name)
                output = os.path.join(self.output_dirs['report_directory'], "multiqc_by_sample", sample.name, f"{sample.name}.RnaSeq.{self.protocol}.multiqc")

                job = multiqc.run(
                        input + "*",
                        output
                        )
                job.name = "multiqc." + sample.name
                job.input_files = self.multiqc_inputs
                job.samples = self.samples
                job.readsets = self.readsets
                jobs.append(job)

        return jobs

    @property
    def step_list(self):
        return self.protocols()[self._protocol]

    def protocols(self):
        return { 'stringtie':
            [
                self.picard_sam_to_fastq,
                self.trimmomatic,
                self.merge_trimmomatic_stats,
                self.sortmerna,
                self.star,
                self.picard_merge_sam_files,
                self.picard_sort_sam,
                self.mark_duplicates,
                self.picard_rna_metrics,
                self.estimate_ribosomal_rna,
                self.rnaseqc2,
                self.wiggle,
                self.raw_counts,
                self.raw_counts_metrics,
                self.stringtie,
                self.stringtie_merge,
                self.stringtie_abund,
                self.ballgown,
                self.differential_expression,
                self.multiqc,
                self.cram_output
            ], 'variants':
            [
                self.picard_sam_to_fastq,
                self.skewer_trimming,
                self.sortmerna,
                self.star,
                self.picard_merge_sam_files,
                self.mark_duplicates,
                self.split_N_trim,
                self.sambamba_merge_splitNtrim_files,
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.recalibration,
                self.gatk_haplotype_caller,
                self.merge_hc_vcf,
                self.run_vcfanno,
                self.variant_filtration,
                self.decompose_and_normalize,
                self.compute_snp_effects,
                self.gemini_annotations,
                self.picard_rna_metrics,
                self.estimate_ribosomal_rna,
                self.rnaseqc2,
                self.gatk_callable_loci,
                self.wiggle,
                self.multiqc,
                self.cram_output
            ], 'cancer':
            [
                self.picard_sam_to_fastq,
                self.skewer_trimming,
                self.sortmerna,
                self.star,
                self.picard_merge_sam_files,
                self.mark_duplicates,
                self.split_N_trim,
                self.sambamba_merge_splitNtrim_files,
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.recalibration,
                self.gatk_haplotype_caller,
                self.merge_hc_vcf,
                self.run_vcfanno,
                self.decompose_and_normalize,
                self.filter_gatk,
                self.report_cpsr,
                self.report_pcgr,
                self.run_star_fusion,
                self.run_arriba,
                self.run_annofuse,
                self.picard_rna_metrics,
                self.estimate_ribosomal_rna,
                self.rnaseqc2,
                self.rseqc,
                self.gatk_callable_loci,
                self.wiggle,
                self.multiqc,
                self.cram_output
            ]
        }


class RnaSeq(RnaSeqRaw):
    __doc__ = RnaSeqRaw.__doc__
    def __init__(self, *args, protocol=None, batch=None, **kwargs):
        self._protocol = protocol
        self.batch_file = batch
        # Add pipeline specific arguments
        super(RnaSeq, self).__init__(*args, **kwargs)

    @classmethod
    def argparser(cls, argparser):
        super().argparser(argparser)
        cls._argparser.add_argument("-t", "--type", help="RNAseq analysis type", dest='protocol',
                                    choices=['stringtie', 'variants', 'cancer'], default="stringtie")
        cls._argparser.add_argument("-b", "--batch", help="batch file (to peform batch effect correction",
                                    type=argparse.FileType('r'))

        return cls._argparser

def main(parsed_args):
    """
    """

    # Pipeline config
    config_files = parsed_args.config

    # Common Pipeline options
    genpipes_file = parsed_args.genpipes_file
    container = parsed_args.container
    clean = parsed_args.clean
    no_json = parsed_args.no_json
    json_pt = parsed_args.json_pt
    force = parsed_args.force
    force_mem_per_cpu = parsed_args.force_mem_per_cpu
    job_scheduler = parsed_args.job_scheduler
    output_dir = parsed_args.output_dir
    steps = parsed_args.steps
    readset_file = parsed_args.readsets_file
    design_file = parsed_args.design_file
    protocol = parsed_args.protocol

    pipeline = RnaSeq(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file, clean=clean, force=force, force_mem_per_cpu=force_mem_per_cpu, job_scheduler=job_scheduler, output_dir=output_dir, design_file=design_file, no_json=no_json, json_pt=json_pt, container=container, protocol=protocol)

    pipeline.submit_jobs()

