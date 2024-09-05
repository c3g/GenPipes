################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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
import logging
import os
import re

# GenPipes Modules
from ...bfx import (
    bash_cmd as bash,
    bcftools,
    bedtools,
    bvatools,
    bwa,
    covseq_tools,
    cutadapt,
    fgbio,
    freebayes,
    gatk4,
    htslib,
    igvtools,
    ivar,
    kraken2,
    multiqc,
    ncovtools,
    picard2,
    qualimap,
    quast,
    sambamba,
    samtools,
    snpeff
    )

from ...core.config import global_conf, _raise, SanitycheckError
from ...core.job import Job, concat_jobs, pipe_jobs
from .. import dnaseq

log = logging.getLogger(__name__)


class CoVSeq(dnaseq.DnaSeqRaw):
    """
CoVSeq Pipeline
================
A pipeline to process and analyze SARS-CoV-2 sequencing data from Illumina platforms. The pipeline uses Cutadapt for adapter trimming, Kraken for taxonomic classification, BWA for read alignment, Sambamba for sorting and indexing, and Freebayes for variant calling. The pipeline also includes a number of metrics to assess the quality.

The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.

Attributes:
    output_dirs (dict): Dictionary of output directories.
    multiqc_inputs (list): List of multiqc input files.
Methods:
    host_reads_removal: Remove host reads from raw reads.
    kraken_analysis: Taxonomic sequence classification system using kraken.
    cutadapt: Raw reads quality trimming and removing adapters.
    mapping_bwa_mem_sambamba: Align filtered reads to a reference genome.
    sambamba_filtering: Filter raw bams with sambamba and an awk cmd to filter by insert size.
    fgbio_trim_primers: Remove primer sequences to individual bam files using fgbio.
    ivar_trim_primers: Remove primer sequences to individual bam files using ivar.
    metrics_sambamba_flagstat: Sambamba flagstsat.
    metrics_bedtools_genomecov: bedtools genome coverage.
    multiple_metrics_picard: Calculates bedtools genome coverage metrics from dnaseq pipeline but on raw AND on filtered bam file.
    metrics_dna_sample_qualimap: Generates metrics with qualimap bamqc.
Parameters:
    protocol (str): Protocol to use for the pipeline.
    """

    def __init__(self, *args, protocol=None, **kwargs):
        if protocol is None:
            self._protocol = 'default'
        # Add pipeline specific arguments
        super(CoVSeq, self).__init__(*args, **kwargs)

    @property
    def output_dirs(self):
        """
        Output directory paths.
        Returns:
            dict: Output directory paths.
        """
        dirs = {
            'raw_reads_directory': os.path.relpath(os.path.join(self.output_dir, 'raw_reads'), self.output_dir),
            'host_removal_directory': os.path.relpath(os.path.join(self.output_dir, "host_removal"), self.output_dir),
            'trim_directory': os.path.relpath(os.path.join(self.output_dir, "cleaned_raw_reads"), self.output_dir),
            'alignment_directory': os.path.relpath(os.path.join(self.output_dir, 'alignment'), self.output_dir),
            'metrics_directory': os.path.relpath(os.path.join(self.output_dir, 'metrics'), self.output_dir),
            'variants_directory': os.path.relpath(os.path.join(self.output_dir, 'variant'), self.output_dir),
            'consensus_directory': os.path.relpath(os.path.join(self.output_dir, 'consensus'), self.output_dir),
            'report_directory': os.path.relpath(os.path.join(self.output_dir, 'report'), self.output_dir)
        }
        return dirs

    @property
    def multiqc_inputs(self):
        """
        List of MultiQC input files.
        Returns:
            list: List of MultiQC input files.
        """
        if not hasattr(self, "_multiqc_inputs"):
            self._multiqc_inputs = []
        return self._multiqc_inputs

    @multiqc_inputs.setter
    def multiqc_inputs(self, value):
        self._multiqc_inputs = value

    def host_reads_removal(self):
        """
        The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
        The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
        BWA output BAM files are then sorted by coordinate using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

        This step takes as input files:
        1. Trimmed FASTQ files if available
        2. Else, FASTQ files from the readset file if available
        3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

        Returns:
            list: List of jobs for host reads removal step.
        """

        jobs = []
        sequencing_center = global_conf.global_get('host_reads_removal', 'sequencing_center', required=False)
        sequencing_technology = global_conf.global_get('host_reads_removal', 'sequencing_technology') if global_conf.global_get('host_reads_removal', 'sequencing_technology', required=False) else "Illumina"
        for readset in self.readsets:
            host_removal_directory = os.path.join(self.output_dirs["host_removal_directory"], readset.sample.name)
            readset_bam = os.path.join(host_removal_directory, f"{readset.name}.hybrid.sorted.bam")
            readset_bam_host_removed_sorted = os.path.join(host_removal_directory, f"{readset.name}.host_removed.sorted.bam")
            readset_bam_host_removed_sorted_index = os.path.join(host_removal_directory, f"{readset.name}.host_removed.sorted.bam.bai")
            readset_bam_host_removed_name_sorted = os.path.join(host_removal_directory, f"{readset.name}.host_removed.nsorted.bam")

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [
                    [readset.fastq1, readset.fastq2]
                ]
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

                output_other = os.path.join(host_removal_directory, f"{readset.name}.host_removed.other.fastq.gz")
                output_single = os.path.join(host_removal_directory, f"{readset.name}.host_removed.single.fastq.gz")
                output_pair1 = os.path.join(host_removal_directory, f"{readset.name}.host_removed.pair1.fastq.gz")
                output_pair2 = os.path.join(host_removal_directory, f"{readset.name}.host_removed.pair2.fastq.gz")

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [
                    [readset.fastq1]
                ]
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                output_pair1 = None
                output_pair2 = None
                output_other = os.path.join(host_removal_directory, f"{readset.name}.host_removed.fastq.gz")
                output_single = None

            else:
                _raise(SanitycheckError(f"""Error: run type "{readset.run_type}" is invalid for readset "{readset.name}" (should be PAIRED_END or SINGLE_END)!"""))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(host_removal_directory),
                        pipe_jobs(
                            [
                                bwa.mem(
                                    fastq1,
                                    fastq2,
                                    read_group="'@RG" + \
                                        f"\\tID:{readset.name}" + \
                                        f"\\tSM:{readset.sample.name}" + \
                                        "\\tLB:" + (readset.library if readset.library else readset.sample.name) + \
                                        ("\\tPU:run" + f"{readset.run}_{readset.lane}" if readset.run and readset.lane else "") + \
                                        (f"\\tCN:{sequencing_center}") + \
                                        (f"\\tPL:{sequencing_technology}") + \
                                        "'",
                                        ini_section='host_reads_removal'
                                ),
                                sambamba.view(
                                    "/dev/stdin",
                                    None,
                                    options="-S -f bam"
                                ),
                                sambamba.sort(
                                    "/dev/stdin",
                                    readset_bam,
                                    tmp_dir=global_conf.global_get('host_reads_removal', 'tmp_dir', required=True),
                                    other_options=global_conf.global_get('host_reads_removal', 'sambamba_sort_other_options', required=False)
                                )
                            ]
                        ),
                        sambamba.view(
                            readset_bam,
                            readset_bam_host_removed_sorted,
                            options=global_conf.global_get('host_reads_removal', 'sambamba_view_other_options')
                        ),
                        sambamba.sort(
                            readset_bam_host_removed_sorted,
                            readset_bam_host_removed_name_sorted,
                            tmp_dir=global_conf.global_get('host_reads_removal', 'tmp_dir', required=True),
                            other_options=global_conf.global_get('host_reads_removal', 'sambamba_name_sort_other_options', required=False)
                        ),
                        sambamba.index(
                            readset_bam_host_removed_sorted,
                            readset_bam_host_removed_sorted_index,
                            other_options=global_conf.global_get('host_reads_removal', 'sambamba_index_other_options', required=False)
                        ),
                        samtools.bam2fq(
                            input_bam=readset_bam_host_removed_name_sorted,
                            output_pair1=output_pair1,
                            output_pair2=output_pair2,
                            output_other=output_other,
                            output_single=output_single,
                            ini_section='host_reads_removal'
                        )
                    ],
                    name=f"host_reads_removal.{readset.name}",
                    samples=[readset.sample]
                )
            )

        return jobs


    def kraken_analysis(self):
        """
        Taxonomic sequence classification system using [kraken](https://github.com/DerrickWood/kraken2).
        Returns:
            list: List of jobs for kraken analysis step.
        """

        # TODO: include kraken analysis and output in metrics
        jobs = []
        for readset in self.readsets:
            host_removal_directory = os.path.join(self.output_dirs["host_removal_directory"], readset.sample.name)
            host_removal_file_prefix = os.path.join(host_removal_directory, readset.name)
            kraken_directory = os.path.join(self.output_dirs["metrics_directory"], "dna", readset.sample.name, "kraken_metrics")
            kraken_out_prefix = os.path.join(kraken_directory, readset.name)

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [
                    [f"{host_removal_file_prefix}.host_removed.pair1.fastq.gz", f"{host_removal_file_prefix}.host_removed.pair2.fastq.gz"]
                ]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append(
                        [readset.fastq1, readset.fastq2]
                    )
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                unclassified_output = [f"{kraken_out_prefix}.unclassified_sequences_1.fastq", f"{kraken_out_prefix}.unclassified_sequences_2.fastq"]
                classified_output = [f"{kraken_out_prefix}.classified_sequences_1.fastq", f"{kraken_out_prefix}.classified_sequences_2.fastq"]

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [
                    [f"{host_removal_file_prefix}.host_removed.single.fastq.gz"]
                ]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                unclassified_output = [f"{kraken_out_prefix}.unclassified_sequences.fastq"]
                classified_output = [f"{kraken_out_prefix}.classified_sequences.fastq"]

            else:
                _raise(SanitycheckError(f"""Error: run type "{readset.run_type}" is invalid for readset "{readset.name}" (should be PAIRED_END or SINGLE_END)!"""))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(kraken_directory),
                        kraken2.kraken2(
                            fastq1,
                            fastq2,
                            kraken_out_prefix,
                            other_options=global_conf.global_get('kraken_analysis', 'kraken2_other_options'),
                            nthread=global_conf.global_get('kraken_analysis', 'kraken2_threads'),
                            database=global_conf.global_get('kraken_analysis', 'kraken2_database')
                        ),
                        bash.pigz(
                            unclassified_output + classified_output,
                            global_conf.global_get('kraken_analysis', 'pigz_threads'),
                            options="-k -f -p",
                            ini_section='kraken_analysis'
                        )
                    ],
                    name=f"kraken_analysis.{readset.name}",
                    samples=[readset.sample],
                    removable_files=unclassified_output + classified_output
                )
            )
        return jobs

    def cutadapt(self):
        """
        Raw reads quality trimming and removing adapters is performed using [Cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html).
        'Adapter1' and 'Adapter2' columns from the readset file ar given to Cutadapt. For PAIRED_END readsets, both adapters are used.
        For SINGLE_END readsets, only Adapter1 is used and left unchanged.
        To trim the front of the read use adapter_5p_fwd and adapter_5p_rev (for PE only) in cutadapt section of ini file.

        This step takes as input files:
        1. FASTQ files from the readset file if available
        2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

        Returns:
            list: List of jobs for cutadapt step.
        """

        jobs = []

        for readset in self.readsets:
            host_removal_directory = os.path.join(self.output_dirs["host_removal_directory"], readset.sample.name)
            host_removal_file_prefix = os.path.join(host_removal_directory, readset.name)
            trim_directory = os.path.join(self.output_dirs['trim_directory'], readset.sample.name)
            trim_file_prefix = os.path.join(trim_directory, readset.name)

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [
                    [f"{host_removal_file_prefix}.host_removed.pair1.fastq.gz", f"{host_removal_file_prefix}.host_removed.pair2.fastq.gz"]
                ]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                adapter_fwd = readset.adapter1
                adapter_rev = readset.adapter2

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [
                    [f"{host_removal_file_prefix}.host_removed.single.fastq.gz"]
                ]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                adapter_fwd = readset.adapter1
                adapter_rev = None

            else:
                _raise(SanitycheckError(f"""Error: run type "{readset.run_type}" is invalid for readset "{readset.name}" (should be PAIRED_END or SINGLE_END)!"""))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(trim_directory),
                        cutadapt.trim(
                            fastq1,
                            fastq2,
                            trim_file_prefix,
                            adapter_fwd,
                            adapter_rev
                        )
                    ],
                    name=f"cutadapt.{readset.name}",
                    samples=[readset.sample]
                )
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

        Returns:
            list: List of jobs for mapping bwa mem sambamba step.
        """

        jobs = []
        sequencing_center = global_conf.global_get('mapping_bwa_mem_sambamba', 'sequencing_center', required=False)
        sequencing_technology = global_conf.global_get('mapping_bwa_mem_sambamba', 'sequencing_technology') if global_conf.global_get('host_reads_removal', 'sequencing_technology', required=False) else "Illumina"
        for readset in self.readsets:
            host_removal_directory = os.path.join(self.output_dirs["host_removal_directory"], readset.sample.name)
            host_removal_file_prefix = os.path.join(host_removal_directory, readset.name)
            trim_directory = os.path.join(self.output_dirs['trim_directory'], readset.sample.name)
            trim_file_prefix = os.path.join(trim_directory, readset.name)
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name)
            readset_bam = os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.bam")
            index_bam = os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.bam.bai")

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [
                    [f"{trim_file_prefix}.trim.pair1.fastq.gz", f"{trim_file_prefix}.trim.pair2.fastq.gz"],
                    [f"{host_removal_file_prefix}.host_removed.pair1.fastq", f"{host_removal_file_prefix}.host_removed.pair2.fastq"]
                ]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [
                    [f"{trim_file_prefix}.trim.single.fastq.gz"],
                    [f"{host_removal_file_prefix}.host_removed.single.fastq"]
                ]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None

            else:
                _raise(SanitycheckError(f"""Error: run type "{readset.run_type}" is invalid for readset "{readset.name}" (should be PAIRED_END or SINGLE_END)!"""))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.dirname(readset_bam)),
                        pipe_jobs(
                            [
                                bwa.mem(
                                    fastq1,
                                    fastq2,
                                    read_group="'@RG" + \
                                        f"\\tID:{readset.name}" + \
                                        f"\\tSM:{readset.sample.name}" + \
                                        "\\tLB:" + (readset.library if readset.library else readset.sample.name) + \
                                        ("\\tPU:run" + f"{readset.run}_{readset.lane}" if readset.run and readset.lane else "") + \
                                        (f"\\tCN:{sequencing_center}") + \
                                        (f"\\tPL:{sequencing_technology}") + \
                                        "'",
                                    ini_section='mapping_bwa_mem_sambamba'
                                ),
                                sambamba.view(
                                    "/dev/stdin",
                                    None,
                                    options=global_conf.global_get('mapping_bwa_mem_sambamba', 'sambamba_view_other_options')
                                ),
                                sambamba.sort(
                                    "/dev/stdin",
                                    readset_bam,
                                    tmp_dir=global_conf.global_get('mapping_bwa_mem_sambamba', 'tmp_dir', required=True),
                                    other_options=global_conf.global_get('mapping_bwa_mem_sambamba', 'sambamba_sort_other_options', required=False)
                                )
                            ]
                        ),
                        sambamba.index(
                            readset_bam,
                            index_bam,
                            other_options=global_conf.global_get('mapping_bwa_mem_sambamba', 'sambamba_index_other_options', required=False)
                        )
                    ],
                    name=f"mapping_bwa_mem_sambamba.{readset.name}",
                    samples=[readset.sample]
                 )
             )

        return jobs

    def sambamba_filtering(self):
        """
        Filter raw bams with [Sambamba](http://lomereiter.github.io/sambamba/index.html) and an awk cmd to filter by insert size
        Returns:
            list: List of jobs for sambamba filtering step.
        """

        library = {}
        for readset in self.readsets:
            # Check the library status
            if not readset.sample in library:
                library[readset.sample] = "SINGLE_END"
            if readset.run_type == "PAIRED_END":
                library[readset.sample] = "PAIRED_END"

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            input_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.bam")
            output_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")

            if library[sample] == "PAIRED_END":
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(os.path.dirname(output_bam)),
                            pipe_jobs(
                                [
                                    sambamba.view(
                                        input_bam=input_bam,
                                        output_bam=None,
                                        options="-h"
                                        ),
                                    Job(
                                        input_files=[],
                                        output_files=[],
                                        command="""awk 'substr($0,1,1)=="@" || \\
  ($9 >= {min_insert_size} && $9 <= {max_insert_size}) || \\
  ($9 <= -{min_insert_size} && $9 >= -{max_insert_size})'""".format(
                                            min_insert_size=global_conf.global_get('sambamba_filtering', 'min_insert_size', required=False),
                                            max_insert_size=global_conf.global_get('sambamba_filtering', 'max_insert_size', required=False)
                                        )
                                    ),
                                    sambamba.view(
                                        input_bam="/dev/stdin",
                                        output_bam=output_bam,
                                        options=global_conf.global_get('sambamba_filtering', 'sambamba_filtering_other_options', required=False)
                                    )
                                ]
                            )
                        ],
                        name=f"sambamba_filtering.{sample.name}",
                        samples=[sample]
                        )
                    )
            else:
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(os.path.dirname(output_bam)),
                            sambamba.view(
                                input_bam=input_bam,
                                output_bam=output_bam,
                                options=global_conf.global_get('sambamba_filtering', 'sambamba_filtering_other_options', required=False)
                            )
                        ],
                        name=f"sambamba_filtering.{sample.name}",
                        samples=[sample]
                    )
                )
        return jobs

    def fgbio_trim_primers(self):
        """
        Remove primer sequences to individual bam files using [fgbio](https://fulcrumgenomics.github.io/fgbio/).
        Returns:
            list: List of jobs for fgbio trim primers step.
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            input_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")
            output_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.primerTrim.bam")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.dirname(output_bam)),
                        fgbio.trim_primers(
                            input_bam,
                            output_bam,
                            hard_clip=True
                        )
                    ],
                    name=f"fgbio_trim_primers.{sample.name}",
                    samples=[sample]
                )
            )
        return jobs

    def ivar_trim_primers(self):
        """
        Remove primer sequences to individual bam files using [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html).
        Returns:
            list: List of jobs for ivar trim primers step.
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            input_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")
            output_prefix = os.path.join(alignment_directory, re.sub(r"\.sorted\.filtered\.bam$", ".primerTrim", os.path.basename(input_bam)))
            output_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.primerTrim.bam")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(alignment_directory),
                        ivar.trim_primers(
                            input_bam,
                            output_prefix
                        ),
                        sambamba.sort(
                            f"{output_prefix}.bam",
                            output_bam,
                            global_conf.global_get('ivar_trim_primers', 'tmp_dir')
                        )
                    ],
                    name=f"ivar_trim_primers.{sample.name}",
                    samples=[sample],
                    removable_files=[f"{output_prefix}.bam"]
                )
            )
        return jobs

    def metrics_sambamba_flagstat(self):
        """
        [Sambamba](http://lomereiter.github.io/sambamba/index.html) flagstsat.
        Returns:
            list: List of jobs for sambamba flagstat step.
        """

        jobs = []
        for sample in self.samples:
            flagstat_directory = os.path.join(self.output_dirs["metrics_directory"], "dna", sample.name, "flagstat")
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            input_bams = [
                os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.primerTrim.bam"),
                os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam"),
                os.path.join(alignment_directory, f"{sample.name}.sorted.bam")
            ]

            for input_bam in input_bams:
                output = os.path.join(flagstat_directory, re.sub(r"\.bam$", ".flagstat", os.path.basename(input_bam)))

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(flagstat_directory),
                            sambamba.flagstat(
                                input_bam,
                                output,
                                global_conf.global_get('sambamba_flagstat', 'flagstat_options')
                            )
                        ],
                        name="sambamba_flagstat." + re.sub(r"\.bam$", "", os.path.basename(input_bam)),
                        samples=[sample]
                    )
                )
        return jobs

    def metrics_bedtools_genomecov(self):
        """
        [bedtools genome coverage](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html).
        Returns:
            list: List of jobs for bedtools genome coverage step.
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            [input_bam] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                ]
            )

            output = os.path.join(alignment_directory, re.sub(r"\.bam$", ".BedGraph", os.path.basename(input_bam)))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(alignment_directory),
                        bedtools.genomecov(
                            input_bam,
                            output
                        )
                    ],
                    name=f"bedtools_genomecov.{sample.name}",
                    samples=[sample]
                )
            )
        return jobs

    def multiple_metrics_picard(self):
        """
        Calculates bedtools genome coverage](https://broadinstitute.github.io/picard/) metrics from dnaseq pipeline but on raw AND on filtered bam file.
        Returns:
            list: List of jobs for multiple metrics picard step.
        """

        # Check the library status
        library = {}
        for readset in self.readsets:
            if not readset.sample in library:
                library[readset.sample] = "SINGLE_END"
            if readset.run_type == "PAIRED_END":
                library[readset.sample] = "PAIRED_END"

        jobs = []
        for sample in self.samples:
            picard_directory = os.path.join(self.output_dirs["metrics_directory"], "dna", sample.name, "picard_metrics")
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            readset = sample.readsets[0]

            input_bams = [
                os.path.join(alignment_directory, f"{sample.name}.sorted.bam"),
                os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")
            ]

            mkdir_job = bash.mkdir(picard_directory, remove=True)

            for input_bam in input_bams:
                jobs.append(
                    concat_jobs(
                        [
                            mkdir_job,
                            gatk4.collect_multiple_metrics(
                                input_bam,
                                os.path.join(picard_directory, re.sub(r"\.bam$", ".all.metrics", os.path.basename(input_bam))),
                                library_type=library[sample]
                            )
                        ],
                        name="multiple_metrics_raw_picard." + re.sub(r"\.bam$", "", os.path.basename(input_bam)),
                        samples=[sample]
                    )
                )
                jobs.append(
                    concat_jobs(
                        [
                            mkdir_job,
                            gatk4.collect_oxog_metrics(
                                input_bam,
                                os.path.join(picard_directory, re.sub(r"\.bam$", ".oxog_metrics.txt", os.path.basename(input_bam)))
                            )
                        ],
                        name="picard_collect_oxog_metrics." + re.sub(r"\.bam$", "", os.path.basename(input_bam)),
                        samples=[sample]
                    )
                )
                jobs.append(
                    concat_jobs(
                        [
                            mkdir_job,
                            gatk4.collect_gcbias_metrics(
                                input_bam,
                                os.path.join(picard_directory, re.sub(r"\.bam$", "", os.path.basename(input_bam)))
                            )
                        ],
                        name="picard_collect_gcbias_metrics." + re.sub(r"\.bam$", "", os.path.basename(input_bam)),
                        samples=[sample]
                    )
                )
        return jobs

    def metrics_dna_sample_qualimap(self):
        """
        Generates metrics with qualimap bamqc:
        BAM QC reports information for the evaluation of the quality of the provided alignment data (a BAM file). 
        In short, the basic statistics of the alignment (number of reads, coverage, GC-content, etc.) are summarized and a number of useful graphs are produced.
        http://qualimap.conesalab.org/doc_html/analysis.html#bamqc
        Returns:
            list: List of jobs for metrics dna sample qualimap step.
        """

        jobs = []
        for sample in self.samples:
            qualimap_directory = os.path.join(self.output_dirs['metrics_directory'], "dna", sample.name, "qualimap", sample.name)
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            [input_file] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.recal.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.matefixed.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.realigned.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                ]
            )
            output = os.path.join(qualimap_directory, "genome_results.txt")

            use_bed = global_conf.global_get('dna_sample_qualimap', 'use_bed', param_type='boolean', required=True)

            options = None
            if use_bed:
                bed = self.samples[0].readsets[0].beds[0]
                options = f"{global_conf.global_get('dna_sample_qualimap', 'qualimap_options')} --feature-file {bed}"

            else:
                options = global_conf.global_get('dna_sample_qualimap', 'qualimap_options')

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            qualimap_directory,
                            remove=False
                        ),
                        bash.mkdir(os.path.join(self.output_dirs['report_directory'], "multiqc_inputs", sample.name)),
                        qualimap.bamqc(
                            input_file,
                            qualimap_directory,
                            output,
                            options,
                            'dna_sample_qualimap'
                        ),
                        bash.ln(
                            os.path.relpath(output, os.path.join(self.output_dirs['report_directory'], "multiqc_inputs", sample.name)),
                            os.path.join(self.output_dirs['report_directory'], "multiqc_inputs", sample.name, os.path.basename(output)),
                            input=output
                        )
                    ],
                    name=f"dna_sample_qualimap.{sample.name}",
                    samples=[sample]
                )
            )
            self.multiqc_inputs.append(output)
        return jobs


    def picard_calculate_hs_metrics(self):
        """
        Compute on target percent of hybridisation based capture.
        Returns:
            list: List of jobs for picard calculate hs metrics step.
        """

        jobs = []

        for sample in self.samples:
            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            if coverage_bed:
                if os.path.isfile(re.sub(r"\.[^.]+$", ".interval_list", coverage_bed)):
                    interval_list = re.sub(r"\.[^.]+$", ".interval_list", coverage_bed)
                else:
                    interval_list = re.sub(r"\.[^.]+$", ".interval_list", os.path.basename(coverage_bed))
                    job = picard2.bed2interval_list(
                        None,
                        coverage_bed,
                        interval_list
                    )
                    job.name = f"interval_list.{os.path.basename(coverage_bed)}"
                    jobs.append(job)

                alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
                [input_file] = self.select_input_files(
                    [
                        [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.recal.bam")],
                        [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.bam")],
                        [os.path.join(alignment_directory, f"{sample.name}.sorted.matefixed.bam")],
                        [os.path.join(alignment_directory, f"{sample.name}.sorted.realigned.bam")],
                        [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                        [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                    ]
                )
                job = gatk4.calculate_hs_metrics(
                    input_file,
                    re.sub("bam$", "onTarget.tsv", input_file),
                    interval_list
                )
                job.name = f"picard_calculate_hs_metrics.{sample.name}"
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
        Returns:
            list: List of jobs for metrics step.
        """

        ##check the library status
        library = {}
        for readset in self.readsets:
            if not readset.sample in library:
                library[readset.sample] = "SINGLE_END"
            if readset.run_type == "PAIRED_END":
                library[readset.sample] = "PAIRED_END"

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            [input_file] = self.select_input_files(
                [
                    # [os.path.join(alignment_directory, f"{sample.name}.sorted.primerTrim.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.recal.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.matefixed.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.realigned.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                ]
            )
            input_file_prefix = re.sub("bam$", "", input_file)

            mkdir_job_normal = bash.mkdir(
                os.path.dirname(input_file_prefix),
                remove=True
            )

            collect_multiple_metrics_normal_job = concat_jobs(
                [
                    mkdir_job_normal,
                    bash.mkdir(os.path.join(self.output_dirs['report_directory'], "multiqc_inputs", sample.name)),
                    gatk4.collect_multiple_metrics(
                        input_file,
                        f"{input_file_prefix}all.metrics",
                        library_type=library[sample]
                    ),
                ]
            )
            for outfile in collect_multiple_metrics_normal_job.report_files:
                self.multiqc_inputs.append(outfile)
                collect_multiple_metrics_normal_job = concat_jobs(
                    [
                        collect_multiple_metrics_normal_job,
                        bash.ln(
                            os.path.relpath(outfile, os.path.join(self.output_dirs['report_directory'], "multiqc_inputs", sample.name)),
                            os.path.join(self.output_dirs['report_directory'], "multiqc_inputs", sample.name, os.path.basename(outfile)),
                            input=outfile
                        )
                    ]
                )
            collect_multiple_metrics_normal_job.name = f"picard_collect_multiple_metrics.{sample.name}"
            collect_multiple_metrics_normal_job.samples = [sample]
            jobs.append(collect_multiple_metrics_normal_job)

            # Compute genome coverage with gatk4
            gatk_depth_of_coverage_job = concat_jobs(
                [
                    mkdir_job_normal,
                    gatk4.depth_of_coverage(
                        input_file,
                        f"{input_file_prefix}all.coverage",
                        bvatools.resolve_readset_coverage_bed(
                            sample.readsets[0]
                        )
                    ),
                ],
                name=f"gatk_depth_of_coverage.{sample.name}.genome",
                samples=[sample]
            )
            jobs.append(gatk_depth_of_coverage_job)

            # Compute genome or target coverage with BVATools
            bvatools_depth_of_coverage_job = concat_jobs(
                [
                    mkdir_job_normal,
                    bvatools.depth_of_coverage(
                        input_file,
                        f"{input_file_prefix}coverage.tsv",
                        bvatools.resolve_readset_coverage_bed(
                            sample.readsets[0]
                        ),
                        other_options=global_conf.global_get('bvatools_depth_of_coverage', 'other_options', required=False)
                    )
                ],
                name=f"bvatools_depth_of_coverage.{sample.name}",
                samples=[sample]
            )
            jobs.append(bvatools_depth_of_coverage_job)

            igvtools_compute_tdf_job = concat_jobs(
                [
                    mkdir_job_normal,
                    igvtools.compute_tdf(
                        input_file,
                        re.sub(r"\.bam$", ".tdf", input_file)
                    )
                ],
                name=f"igvtools_compute_tdf.{sample.name}",
                samples=[sample]
            )
            jobs.append(igvtools_compute_tdf_job)

        return jobs

    def covseq_metrics(self):
        """
        Gathering multiple metrics.
        Returns:
            list: List of jobs for multiple metrics step.
        """

        jobs = []

        jobs.extend(self.multiple_metrics_picard())
        jobs.extend(self.metrics_dna_sample_qualimap())
        jobs.extend(self.metrics_sambamba_flagstat())
        jobs.extend(self.metrics_bedtools_genomecov())
        jobs.extend(self.picard_calculate_hs_metrics())
        jobs.extend(self.metrics())

        for job in jobs:
            self.multiqc_inputs.extend(job.output_files)

        return jobs

    def ivar_calling(self):

        """
        [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) calling variants.
        Returns:
            list: List of jobs for ivar calling step.
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            variant_directory = os.path.join(self.output_dirs["variants_directory"], sample.name)

            [input_bam] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.primerTrim.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                ]
            )

            output_prefix = f"{os.path.join(variant_directory, sample.name)}.variants"
            output_tsv = f"{output_prefix}.tsv"
            output_vcf = os.path.join(variant_directory, re.sub(r"\.bam$", "", os.path.basename(input_bam))) + ".vcf"

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(variant_directory),
                        pipe_jobs(
                            [
                                samtools.mpileup(
                                    input_bam,
                                    output=None,
                                    region=None,
                                    regionFile=None,
                                    ini_section='ivar_call_variants'
                                ),
                                ivar.call_variants(output_prefix)
                            ]
                        ),
                        ivar.tsv_to_vcf(
                            output_tsv,
                            output_vcf
                        ),
                        htslib.bgzip_tabix(
                            output_vcf,
                            f"{output_vcf}.gz"
                        )
                    ],
                    name=f"ivar_call_variants.{sample.name}",
                    samples=[sample],
                    removable_files=[output_vcf]
                )
            )
        return jobs

    def freebayes_calling(self):
        """
        [FreeBayes](https://github.com/freebayes/freebayes) is a haplotype-based variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.

        This method avoids one of the core problems with alignment-based variant detection that identical sequences may have multiple possible alignments.
        Returns:
            list: List of jobs for freebayes calling step.
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            variant_directory = os.path.join(self.output_dirs["variants_directory"], sample.name)

            [input_bam] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.primerTrim.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                ]
            )

            output_prefix = os.path.join(variant_directory, f"{sample.name}.freebayes_calling")
            output_gvcf = f"{output_prefix}.gvcf"

            output_masks = f"{output_prefix}.mask.txt"
            output_variants = f"{output_prefix}.variants.vcf"
            output_consensus = f"{output_prefix}.consensus.vcf"
            output_consensus_norm = re.sub(r"\.vcf$", ".norm.vcf.gz", output_consensus)

            job = concat_jobs(
                [
                    bash.mkdir(variant_directory),
                    freebayes.freebayes(
                        input_bam,
                        output_gvcf,
                        options="--gvcf --gvcf-dont-use-chunk true",
                        ini_section='freebayes_call_variants'
                    ),
                    freebayes.process_gvcf(
                        output_gvcf,
                        output_masks,
                        output_variants,
                        output_consensus,
                        ini_section='freebayes_call_variants'
                    )
                ]
            )
            for file in [output_variants, output_consensus]:
                output_vcf_gz = re.sub(r"\.vcf$", ".norm.vcf.gz", file)
                job = concat_jobs(
                    [
                        job,
                        pipe_jobs(
                            [
                                bcftools.norm(
                                    file,
                                    None,
                                    "-f " + global_conf.global_get("DEFAULT", 'genome_fasta', param_type='filepath'),
                                    ini_section='freebayes_call_variants'
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    output_vcf_gz
                                )
                            ]
                        )
                    ]
                )
            for flag in ["ambiguous", "fixed"]:
                output_vcf_gz = re.sub(r"\.consensus\.norm\.vcf\.gz$", f".{flag}.norm.vcf.gz", output_consensus_norm)
                job = concat_jobs(
                    [
                        job,
                        pipe_jobs(
                            [
                                bash.cat(
                                    output_consensus_norm,
                                    None,
                                    zip=True
                                ),
                                bash.awk(
                                    None,
                                    None,
                                    f"-v vartag=ConsensusTag={flag} '$0 ~ /^#/ || $0 ~ vartag'"
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    output_vcf_gz
                                )
                            ]
                        )
                    ]
                )
            job.name = f"freebayes_call_variants.{sample.name}"
            job.samples = [sample]
            job.removable_files = [output_gvcf]
            jobs.append(job)

        return jobs

    def snpeff_annotate(self):
        """
        Consensus annotation with [SnpEff](https://pcingola.github.io/SnpEff/).
        Returns:
            list: List of jobs for snpeff annotate step.
        """

        jobs = []

        for sample in self.samples:
            variant_directory = os.path.join(self.output_dirs["variants_directory"], sample.name)
            metrics_ivar_prefix = os.path.join(self.output_dirs["metrics_directory"], "dna", sample.name, "snpeff_metrics_ivar", f"{sample.name}.snpEff")
            metrics_freebayes_prefix = os.path.join(self.output_dirs["metrics_directory"], "dna", sample.name, "snpeff_metrics_freebayes", f"{sample.name}.snpEff")

            [ivar_vcf] = self.select_input_files(
                [
                    [os.path.join(variant_directory, f"{sample.name}.sorted.filtered.primerTrim.vcf.gz")],
                    [os.path.join(variant_directory, f"{sample.name}.sorted.filtered.vcf.gz")],
                    [os.path.join(variant_directory, f"{sample.name}.sorted.vcf.gz")]
                ]
            )

            ivar_output_vcf = os.path.join(variant_directory, re.sub(r"\.vcf\.gz$", ".annotate.vcf", os.path.basename(ivar_vcf)))

            freebayes_vcf = os.path.join(variant_directory, f"{sample.name}.freebayes_calling.fixed.norm.vcf.gz")
            freebayes_output_vcf = os.path.join(variant_directory, re.sub(r"\.vcf\.gz$", ".annotate.vcf", os.path.basename(freebayes_vcf)))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.dirname(metrics_ivar_prefix)),
                        snpeff.snpeff_annotate(
                            ivar_vcf,
                            ivar_output_vcf,
                            metrics_ivar_prefix
                        ),
                        htslib.bgzip_tabix(
                            ivar_output_vcf,
                            f"{ivar_output_vcf}.gz"
                        ),
                        bash.mkdir(os.path.dirname(metrics_freebayes_prefix)),
                        snpeff.snpeff_annotate(
                            freebayes_vcf,
                            freebayes_output_vcf,
                            metrics_freebayes_prefix
                        ),
                        htslib.bgzip_tabix(
                            freebayes_output_vcf,
                            f"{freebayes_output_vcf}.gz"
                        )
                    ],
                    name=f"snpeff_annotate.{sample.name}",
                    samples=[sample]
                )
            )
        return jobs

    def ivar_create_consensus(self):
        """
        Create consensus with [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) through a [samtools](http://samtools.sourceforge.net/) mpileup
        Returns:
            list: List of jobs for ivar create consensus step.
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            consensus_directory = os.path.join(self.output_dirs["consensus_directory"], sample.name)
            [input_bam] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.primerTrim.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                ]
            )

            output_prefix = os.path.join(consensus_directory, re.sub(r"\.bam$", ".consensus", os.path.basename(input_bam)))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.dirname(output_prefix)),
                        pipe_jobs(
                            [
                                samtools.mpileup(
                                    input_bam,
                                    output=None,
                                    region=None,
                                    regionFile=None,
                                    ini_section='ivar_create_consensus'
                                ),
                                ivar.create_consensus(
                                    output_prefix
                                )
                            ]
                        )
                    ],
                    name=f"ivar_create_consensus.{sample.name}",
                    samples=[sample],
                    removable_files=[f"{output_prefix}.fa"]
                )
            )
        return jobs

    def bcftools_create_consensus(self):
        """
        [bcftools](https://samtools.github.io/bcftools/bcftools.html) consensus creation
        Returns:
            list: List of jobs for bcftools create consensus step.
        """

        jobs = []
        for sample in self.samples:
            consensus_directory = os.path.join(self.output_dirs["consensus_directory"], sample.name)
            variant_directory = os.path.join(self.output_dirs["variants_directory"], sample.name)

            input_prefix = os.path.join(variant_directory, f"{sample.name}.freebayes_calling")

            intput_masks = f"{input_prefix}.mask.txt"
            input_ambiguous_norm = f"{input_prefix}.ambiguous.norm.vcf.gz"
            input_fixed_norm = f"{input_prefix}.fixed.norm.vcf.gz"

            output_prefix = os.path.join(consensus_directory, f"{sample.name}.freebayes_calling")
            output_ambiguous_fasta = f"{output_prefix}.ambiguous.fasta"
            output_consensus_fasta = f"{output_prefix}.consensus.fasta"

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.dirname(output_prefix)),
                        bcftools.consensus(
                            input_ambiguous_norm,
                            output_ambiguous_fasta,
                            f"-f {global_conf.global_get('DEFAULT', 'genome_fasta', param_type='filepath')} -I "
                        ),
                        pipe_jobs(
                            [
                                bcftools.consensus(
                                    input_fixed_norm,
                                    None,
                                    f"-f {output_ambiguous_fasta} -m {intput_masks}"
                                    ),
                                bash.sed(
                                    None,
                                    output_consensus_fasta,
                                    f"""s/{global_conf.global_get("DEFAULT", 'assembly_synonyms')}/{sample.name}/"""
                                )
                            ]
                        )
                    ],
                    name=f"bcftools_create_consensus.{sample.name}",
                    samples=[sample]
                )
            )
        return jobs

    def quast_consensus_metrics(self):
        """
        Generate [QUAST](http://quast.sourceforge.net/) metrics on consensus
        Returns:
            list: List of jobs for quast consensus metrics step.
        """

        jobs = []
        for sample in self.samples:
            consensus_directory = os.path.join(self.output_dirs["consensus_directory"], sample.name)
            ivar_output_dir = os.path.join(self.output_dirs["metrics_directory"], "dna", sample.name, "quast_metrics_ivar")
            [ivar_consensus] = self.select_input_files(
                [
                    [os.path.join(consensus_directory, f"{sample.name}.sorted.filtered.primerTrim.consensus.fa")],
                    [os.path.join(consensus_directory, f"{sample.name}.sorted.filtered.consensus.fa")],
                    [os.path.join(consensus_directory, f"{sample.name}.sorted.consensus.fa")]
                ]
            )
            freebayes_output_dir = os.path.join(self.output_dirs["metrics_directory"], "dna", sample.name, "quast_metrics_freebayes")
            freebayes_consensus = os.path.join(consensus_directory, f"{sample.name}.freebayes_calling.consensus.fasta")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(ivar_output_dir),
                        quast.quast(
                            ivar_consensus,
                            ivar_output_dir
                        ),
                        bash.mkdir(freebayes_output_dir),
                        quast.quast(
                            freebayes_consensus,
                            freebayes_output_dir
                        )
                    ],
                    name=f"quast_consensus_metrics.{sample.name}",
                    samples=[sample]
                )
            )
        return jobs

    def rename_consensus_header_ivar(self):
        """
        Rename reads headers after [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) variants calling.
        Returns:
            list: List of jobs for rename consensus header ivar step.
        """

        jobs = []
        country = global_conf.global_get('rename_consensus_header', 'country', required=False)
        province = global_conf.global_get('rename_consensus_header', 'province', required=False)
        year = global_conf.global_get('rename_consensus_header', 'year', required=False)
        seq_method = global_conf.global_get('rename_consensus_header', 'seq_method', required=False)

        for sample in self.samples:
            consensus_directory = os.path.join(self.output_dirs["consensus_directory"], sample.name)
            [ivar_consensus] = self.select_input_files(
                [
                    [os.path.join(consensus_directory, f"{sample.name}.sorted.filtered.primerTrim.consensus.fa")],
                    [os.path.join(consensus_directory, f"{sample.name}.sorted.filtered.consensus.fa")],
                    [os.path.join(consensus_directory, f"{sample.name}.sorted.consensus.fa")]
                ]
            )
            quast_ivar_directory = os.path.join(self.output_dirs["metrics_directory"], "dna", sample.name, "quast_metrics_ivar")
            quast_ivar_html = os.path.join(quast_ivar_directory, "report.html")
            quast_ivar_tsv = os.path.join(quast_ivar_directory, "report.tsv")

            ivar_output_fa = os.path.join(consensus_directory, f"{sample.name}.consensus.fasta")
            ivar_output_status_fa = os.path.join(consensus_directory, f"{sample.name}.consensus.{global_conf.global_get('rename_consensus_header', 'sequencing_technology', required=True)}.${{IVAR_STATUS}}.fasta")

            variant_directory = os.path.join(self.output_dirs["variants_directory"], sample.name)
            [ivar_vcf] = self.select_input_files(
                [
                    [os.path.join(variant_directory, f"{sample.name}.sorted.filtered.primerTrim.vcf.gz")],
                    [os.path.join(variant_directory, f"{sample.name}.sorted.filtered.vcf.gz")],
                    [os.path.join(variant_directory, f"{sample.name}.sorted.vcf.gz")]
                ]
            )
            ivar_annotated_vcf = os.path.join(variant_directory, re.sub(r"\.vcf\.gz$", ".annotate.vcf", os.path.basename(ivar_vcf)))

            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            [input_bam] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                ]
            )

            bedgraph_file = os.path.join(alignment_directory, re.sub(r"\.bam$", ".BedGraph", os.path.basename(input_bam)))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.dirname(ivar_output_fa)),
                        Job(
                            input_files=[quast_ivar_tsv, quast_ivar_html],
                            output_files=[],
                            command=f"""\\
ivar_cons_len=`grep -oP "Total length \\(>= 0 bp\\)\\t\\K.*?(?=$)" {quast_ivar_tsv}`
ivar_N_count=`grep -oP "# N's\\",\\"quality\\":\\"Less is better\\",\\"values\\":\\[\\K.*?(?=])" {quast_ivar_html}`
ivar_cons_perc_N=`echo "scale=2; 100*$ivar_N_count/$ivar_cons_len" | bc -l`
ivar_frameshift=`if grep -q "frameshift_variant" {ivar_annotated_vcf}; then echo "FLAG"; fi`
genome_size=`awk '{{print $2}}' {global_conf.global_get('DEFAULT', 'igv_genome', required=True)}`
bam_cov50X=`awk '{{if ($4 > 50) {{count = count + $3-$2}}}} END {{if (count) {{print count}} else {{print 0}}}}' {bedgraph_file}`
bam_cov50X=`echo "scale=2; 100*$bam_cov50X/$genome_size" | bc -l`
IVAR_STATUS=`awk -v bam_cov50X=$bam_cov50X -v ivar_frameshift=$ivar_frameshift -v ivar_cons_perc_N=$ivar_cons_perc_N 'BEGIN {{ if (ivar_cons_perc_N < 5 && ivar_frameshift != "FLAG" && bam_cov50X >= 90) {{print "pass"}}  else if (ivar_cons_perc_N > 10) {{print "rej"}} else if ((ivar_cons_perc_N >= 5 && ivar_cons_perc_N <= 10) || ivar_frameshift == "FLAG" || bam_cov50X < 90) {{print "flag"}} }}'`
export IVAR_STATUS"""
                        ),
                        Job(
                            input_files=[ivar_consensus],
                            output_files=[ivar_output_fa],
                            command=f"""\\
awk '/^>/{{print ">{country}/{province}-{sample.name}/{year} seq_method:{seq_method}|assemb_method:ivar|snv_call_method:ivar"; next}}{{print}}' < {ivar_consensus} > {ivar_output_status_fa} && \\
ln -sf {os.path.basename(ivar_output_status_fa)} {ivar_output_fa}"""
                        )
                    ],
                    name=f"rename_consensus_header.{sample.name}",
                    samples=[sample]
                )
            )
        return jobs

    def rename_consensus_header_freebayes(self):
        """
        Rename reads headers after [FreeBayes](https://github.com/freebayes/freebayes) variants calling.
        Returns:
            list: List of jobs for rename consensus header freebayes step.
        """

        jobs = []
        country = global_conf.global_get('rename_consensus_header', 'country', required=False)
        province = global_conf.global_get('rename_consensus_header', 'province', required=False)
        year = global_conf.global_get('rename_consensus_header', 'year', required=False)
        seq_method = global_conf.global_get('rename_consensus_header', 'seq_method', required=False)

        for sample in self.samples:
            consensus_directory = os.path.join(self.output_dirs["consensus_directory"], sample.name)
            freebayes_consensus = os.path.join(consensus_directory, f"{sample.name}.freebayes_calling.consensus.fasta")

            quast_freebayes_directory = os.path.join(self.output_dirs["metrics_directory"], "dna", sample.name, "quast_metrics_freebayes")
            quast_freebayes_html = os.path.join(quast_freebayes_directory, "report.html")
            quast_freebayes_tsv = os.path.join(quast_freebayes_directory, "report.tsv")

            freebayes_output_fa = os.path.join(consensus_directory, f"{sample.name}.freebayes_calling.consensus.renamed.fasta")
            freebayes_output_status_fa = os.path.join(consensus_directory, f"{sample.name}.freebayes_calling.consensus.{global_conf.global_get('rename_consensus_header', 'sequencing_technology', required=True)}.${{FREEBAYES_STATUS}}.fasta")

            variant_directory = os.path.join(self.output_dirs["variants_directory"], sample.name)
            freebayes_vcf = os.path.join(variant_directory, f"{sample.name}.freebayes_calling.fixed.norm.vcf.gz")
            freebayes_annotated_vcf = os.path.join(variant_directory, re.sub(r"\.vcf\.gz$", ".annotate.vcf", os.path.basename(freebayes_vcf)))

            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            [input_bam] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                ]
            )

            bedgraph_file = os.path.join(alignment_directory, re.sub(r"\.bam$", ".BedGraph", os.path.basename(input_bam)))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.dirname(freebayes_output_fa)),
                        Job(
                            input_files=[quast_freebayes_tsv, quast_freebayes_html],
                            output_files=[],
                            command=f"""\\
freebayes_cons_len=`grep -oP "Total length \\(>= 0 bp\\)\\t\\K.*?(?=$)" {quast_freebayes_tsv}`
freebayes_N_count=`grep -oP "# N's\\",\\"quality\\":\\"Less is better\\",\\"values\\":\\[\\K.*?(?=])" {quast_freebayes_html}`
freebayes_cons_perc_N=`echo "scale=2; 100*$freebayes_N_count/$freebayes_cons_len" | bc -l`
freebayes_frameshift=`if grep -q "frameshift_variant" {freebayes_annotated_vcf}; then echo "FLAG"; fi`
genome_size=`awk '{{print $2}}' {global_conf.global_get('DEFAULT', 'igv_genome', required=False)}`
bam_cov50X=`awk '{{if ($4 > 50) {{count = count + $3-$2}}}} END {{if (count) {{print count}} else {{print 0}}}}' {bedgraph_file}`
bam_cov50X=`echo "scale=2; 100*$bam_cov50X/$genome_size" | bc -l`
FREEBAYES_STATUS=`awk -v bam_cov50X=$bam_cov50X -v freebayes_frameshift=$freebayes_frameshift -v freebayes_cons_perc_N=$freebayes_cons_perc_N 'BEGIN {{ if (freebayes_cons_perc_N < 5 && freebayes_frameshift != "FLAG" && bam_cov50X >= 90) {{print "pass"}}  else if (freebayes_cons_perc_N > 10) {{print "rej"}} else if ((freebayes_cons_perc_N >= 5 && freebayes_cons_perc_N <= 10) || freebayes_frameshift == "FLAG" || bam_cov50X < 90) {{print "flag"}} }}'`
export FREEBAYES_STATUS"""
                        ),
                        Job(
                            input_files=[freebayes_consensus],
                            output_files=[freebayes_output_fa],
                            command=f"""\\
awk '/^>/{{print ">{country}/{province}-{sample.name}/{year} seq_method:{seq_method}|assemb_method:bcftools|snv_call_method:freebayes"; next}}{{print}}' < {freebayes_consensus} > {freebayes_output_status_fa} && \\
ln -sf {os.path.basename(freebayes_output_status_fa)} {freebayes_output_fa}"""
                        )
                    ],
                    name=f"rename_consensus_header.{sample.name}",
                    samples=[sample]
                )
            )
        return jobs

    def run_multiqc(self):
        """
        Run [multiqc](https://multiqc.info/) on all samples.
        Returns:
            list: List of jobs for multiqc step.
        """

        jobs = []

        metrics_directory = os.path.join(self.output_dirs["metrics_directory"], "dna")
        output = os.path.join(metrics_directory, "CoVSeq.multiqc_report")

        job = multiqc.run(
            self.multiqc_inputs,
            output
        )
        job.name = "multiqc_all_samples"
        job.samples = self.samples

        jobs.append(job)

        return jobs

    def ncovtools_quickalign(self):
        """
        Uses [ncov-tools](https://github.com/jts/ncov-tools) quickalign to provides summary statistics, which can be used to determine the sequencing quality and evolutionary novelty of input genomes (e.g. number of new mutations and indels). 

        It uses ivar consensus as well as freebayes consensus to arrive at the alignment decisions.
        Returns:
            list: List of jobs for ncovtools quickalign step.
        """

        jobs = []

        for sample in self.samples:
            ncovtools_quickalign_directory = os.path.join(self.output_dirs["metrics_directory"], "dna", sample.name, "ncovtools_quickalign")
            output = os.path.join(ncovtools_quickalign_directory, f"{sample.name}_ivar_vs_freebayes.vcf")
            consensus_directory = os.path.join(self.output_dirs["consensus_directory"], sample.name)
            ivar_consensus = os.path.join(consensus_directory, f"{sample.name}.consensus.fasta")
            freebayes_consensus = os.path.join(consensus_directory, f"{sample.name}.freebayes_calling.consensus.renamed.fasta")
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(ncovtools_quickalign_directory),
                        Job(
                            input_files=[ivar_consensus, freebayes_consensus],
                            output_files=[output],
                            module_entries=[
                                ['ncovtools_quickalign', 'module_python'],
                                ['ncovtools_quickalign', 'module_ncov_random_scripts'],
                            ],
                            command=f"""\\
quick_align.py -r {ivar_consensus} -g {freebayes_consensus} -o vcf > {output}"""
                        )
                    ],
                    name=f"ncovtools_quickalign.{sample.name}",
                    samples=[sample]
                )
            )
        return jobs

    def prepare_table(self):
        """
        Gathers all analysis data for [QUAST](http://quast.sourceforge.net/), [kraken](https://github.com/DerrickWood/kraken2) and other metrics and module details.
        Returns:
            list: List of jobs for prepare table step.
        """
        jobs = []

        metrics_directory = os.path.join(self.output_dirs["metrics_directory"], "dna")

        readset_file=os.path.relpath(self.readsets_file.name, self.output_dir)

        run_metadata = os.path.join(self.output_dirs["report_directory"], "run_metadata.csv")
        software_version = os.path.join(self.output_dirs["report_directory"], "software_versions.csv")

        modules = []
        # Retrieve all unique module version values in config files
        # assuming that all module key names start with "module_"
        for section in global_conf.sections():
            for name, value in global_conf.items(section):
                if re.search("^module_", name) and value not in modules:
                    modules.append(value)

        # Finding all kraken outputs
        kraken_outputs = []
        library = {}
        for readset in self.readsets:
            ##check the library status
            if not readset.sample in library:
                library[readset.sample] = "SINGLE_END"
            if readset.run_type == "PAIRED_END":
                library[readset.sample] = "PAIRED_END"

            kraken_directory = os.path.join(self.output_dirs["metrics_directory"], "dna", readset.sample.name, "kraken_metrics")
            kraken_out_prefix = os.path.join(kraken_directory, readset.name)
            kraken_output = f"{kraken_out_prefix}.kraken2_output"
            kraken_report = f"{kraken_out_prefix}.kraken2_report"
            kraken_outputs.extend((kraken_output, kraken_report))
            if readset.run_type == "PAIRED_END":
                unclassified_output_1 = f"{kraken_out_prefix}.unclassified_sequences_1.fastq"
                unclassified_output_2 = f"{kraken_out_prefix}.unclassified_sequences_2.fastq"
                classified_output_1 = f"{kraken_out_prefix}.classified_sequences_1.fastq"
                classified_output_2 = f"{kraken_out_prefix}.classified_sequences_2.fastq"
                kraken_outputs.extend((unclassified_output_1, unclassified_output_2, classified_output_1, classified_output_2))

            elif readset.run_type == "SINGLE_END":
                unclassified_output = f"{kraken_out_prefix}.unclassified_sequences.fastq"
                classified_output = f"{kraken_out_prefix}.classified_sequences.fastq"
                kraken_outputs.extend((unclassified_output, classified_output))

        quast_outputs = []
        flagstat_outputs = []
        bedgraph_outputs = []
        picard_outputs = []
        run_name = global_conf.global_get('prepare_report', 'run_name', required=True)
        cluster_server = global_conf.global_get('prepare_report', 'cluster_server')
        assembly_synonyms = global_conf.global_get('prepare_report', 'assembly_synonyms')
        sequencing_technology = global_conf.global_get('prepare_report', 'sequencing_technology')
        for sample in self.samples:
            # Finding quast outputs
            ivar_quast_output_dir = os.path.join(self.output_dirs["metrics_directory"], "dna", sample.name, "quast_metrics_ivar")
            quast_outputs.extend(
                [
                    os.path.join(ivar_quast_output_dir, "report.html"),
                    os.path.join(ivar_quast_output_dir, "report.pdf"),
                    os.path.join(ivar_quast_output_dir, "report.tex"),
                    os.path.join(ivar_quast_output_dir, "report.tsv"),
                    os.path.join(ivar_quast_output_dir, "report.txt")
                ]
            )
            # Finding flagstat outputs
            flagstat_directory = os.path.join(self.output_dirs["metrics_directory"], "dna", sample.name, "flagstat")
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            input_bams = [
                os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.primerTrim.bam"),
                os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam"),
                os.path.join(alignment_directory, f"{sample.name}.sorted.bam")
            ]
            for input_bam in input_bams:
                flagstat_outputs.append(os.path.join(flagstat_directory, re.sub(r"\.bam$", ".flagstat", os.path.basename(input_bam))))
            # Finding BedGraph file
            [input_bam] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                ]
            )
            bedgraph_outputs.append(os.path.join(alignment_directory, re.sub(r"\.bam$", ".BedGraph", os.path.basename(input_bam))))

            # Picard collect multiple metrics outputs
            picard_directory = os.path.join(self.output_dirs["metrics_directory"], "dna", sample.name, "picard_metrics")
            picard_out = os.path.join(picard_directory, re.sub(r"\.bam$", "", os.path.basename(input_bam)) + ".all.metrics")

            if library[sample] == "PAIRED_END":
                picard_outputs.extend(
                    [
                        f"{picard_out}.alignment_summary_metrics",
                        f"{picard_out}.insert_size_metrics"
                    ]
                )
            else:
                picard_outputs.extend(
                    [
                        f"{picard_out}.alignment_summary_metrics",
                    ]
                )

        covid_collect_metrics_inputs = []
        covid_collect_metrics_inputs.extend(kraken_outputs)
        covid_collect_metrics_inputs.append(input_bam)
        covid_collect_metrics_inputs.extend(quast_outputs)
        covid_collect_metrics_inputs.extend(flagstat_outputs)
        covid_collect_metrics_inputs.extend(bedgraph_outputs)
        covid_collect_metrics_inputs.extend(picard_outputs)
        modules_list = "\n".join(modules)
        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(metrics_directory),
                    bash.mkdir(os.path.dirname(run_metadata)),
                    covseq_tools.covid_collect_metrics(
                        readset_file,
                        covid_collect_metrics_inputs
                    ),
                    Job(
                        input_files=[],
                        output_files=[run_metadata, software_version],
                        module_entries=[],
                        command=f"""\\
echo "Preparing to run metadata..." && \\
echo "run_name,{run_name}
genpipes_version,{self.genpipes_version}
cluster_server,{cluster_server}
assembly_synonyms,{assembly_synonyms}
sequencing_technology,{sequencing_technology}" > {run_metadata} && \\
echo "Software Versions
{modules_list}" > {software_version}"""
                    )
                ],
                name="prepare_table." + global_conf.global_get('prepare_report', 'run_name', required=True)
            )
        )
        return jobs

    def prepare_report_ivar(self):
        """
        Prepare [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) analysis report.
        Returns:
            list: List of jobs for prepare report ivar step.
        """
        jobs = []

        readset_file = os.path.relpath(self.readsets_file.name, self.output_dir)
        ivar_readset_file_report = os.path.join(self.output_dirs["report_directory"], "report.readset_ivar.tsv")

        software_version = os.path.join(self.output_dirs["report_directory"], "software_versions.csv")
        run_metadata = os.path.join(self.output_dirs["report_directory"], "run_metadata.csv")

        ivar_ncovtools_directory = os.path.join(self.output_dirs["report_directory"], "ncov_tools_ivar")
        ivar_metadata = os.path.join(ivar_ncovtools_directory, "metadata.tsv")
        ivar_ncovtools_data_directory = os.path.join(ivar_ncovtools_directory, "data")
        ivar_ncovtools_config = os.path.join(ivar_ncovtools_directory, "config.yaml")

        job = concat_jobs(
            [
                bash.mkdir(ivar_ncovtools_data_directory),
                Job(
                    input_files=[],
                    output_files=[
                        ivar_readset_file_report,
                        ivar_metadata
                    ],
                    command=f"""\\
head -n 1 {readset_file} > {ivar_readset_file_report} && \\
echo -e "sample\\tct\\tdate" > {ivar_metadata}"""
                )
            ]
        )

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            filtered_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")
            primer_trimmed_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.primerTrim.bam")
            consensus_directory = os.path.join(self.output_dirs["consensus_directory"], sample.name)
            ivar_consensus = os.path.join(consensus_directory, f"{sample.name}.consensus.fasta")
            variant_directory = os.path.join(self.output_dirs["variants_directory"], sample.name)
            ivar_variants = os.path.join(variant_directory, f"{sample.name}.variants.tsv")

            ivar_output_filtered_bam = os.path.join(ivar_ncovtools_data_directory, os.path.basename(filtered_bam))
            ivar_output_primer_trimmed_bam = os.path.join(ivar_ncovtools_data_directory, os.path.basename(primer_trimmed_bam))
            ivar_output_consensus = os.path.join(ivar_ncovtools_data_directory, os.path.basename(ivar_consensus))
            ivar_output_variants = os.path.join(ivar_ncovtools_data_directory, os.path.basename(ivar_variants))

            job = concat_jobs(
                [
                    job,
                    Job(
                        input_files=[
                            filtered_bam,
                            primer_trimmed_bam,
                            ivar_consensus,
                            ivar_variants
                        ],
                        output_files=[
                            ivar_output_filtered_bam,
                            ivar_output_primer_trimmed_bam,
                            ivar_output_consensus,
                            ivar_output_variants
                        ],
                        command=f"""\\
echo "Linking files for ncov_tools for sample {sample.name}..." && \\
if [ "$(ls -1 {filtered_bam})" != "" ] && [ "$(ls -1 {primer_trimmed_bam})" != "" ] && [ "$(ls -1 {ivar_consensus})" != "" ] && [ "$(ls -1 {ivar_variants})" != "" ];
  then
    ln -fs {os.path.relpath(filtered_bam, os.path.dirname(ivar_output_filtered_bam))} {ivar_output_filtered_bam} && \\
    ln -fs {os.path.relpath(primer_trimmed_bam, os.path.dirname(ivar_output_primer_trimmed_bam))} {ivar_output_primer_trimmed_bam} && \\
    ln -fs {os.path.relpath(ivar_consensus, os.path.dirname(ivar_output_consensus))} {ivar_output_consensus} && \\
    ln -fs {os.path.relpath(ivar_variants, os.path.dirname(ivar_output_variants))} {ivar_output_variants} && \\
    grep -w {sample.name} {readset_file} >> {ivar_readset_file_report} && \\
    echo -e "{sample.name}\\tNA\\tNA" >> {ivar_metadata}
fi"""
                    )
                ],
                samples=[sample]
            )
        jobs.append(
            concat_jobs(
                [
                    job,
                    ncovtools.run_ncovtools(
                        ivar_output_filtered_bam,
                        ivar_output_primer_trimmed_bam,
                        ivar_output_consensus,
                        ivar_output_variants,
                        readset_file,
                        ivar_metadata,
                        ivar_ncovtools_directory,
                        ivar_ncovtools_config,
                        self.output_dir
                    ),
                    Job(
                        input_files=[],
                        output_files=[],
                        module_entries=[
                            ['prepare_report', 'module_R'],
                            ['prepare_report', 'module_CoVSeQ_tools'],
                            ['prepare_report', 'module_pandoc']
                        ],
                        command=f"""\
module purge && \\
module load {global_conf.global_get('prepare_report', 'module_R')} {global_conf.global_get('prepare_report', 'module_CoVSeQ_tools')} {global_conf.global_get('prepare_report', 'module_pandoc')}"""
                    ),
                    covseq_tools.generate_report_tables(
                        ivar_readset_file_report,
                        output_name_pattern=os.path.join("report", "report_metrics_ivar")
                    ),
                    covseq_tools.render_report(
                        software_version,
                        run_metadata,
                        output_name_pattern=os.path.join("report", "report_metrics_ivar"),
                        caller="ivar"
                    )
                ],
                name=f"prepare_report.{global_conf.global_get('prepare_report', 'run_name', required=True)}"
            )
        )

        return jobs

    def prepare_report_freebayes(self):
        """
        Prepare [FreeBayes](https://github.com/freebayes/freebayes) analysis report.
        Returns:
            list: List of jobs for prepare report freebayes step.
        """

        jobs = []

        readset_file=os.path.relpath(self.readsets_file.name, self.output_dir)
        freebayes_readset_file_report=os.path.join(self.output_dirs["report_directory"], "report.readset_freebayes.tsv")

        software_version = os.path.join(self.output_dirs["report_directory"], "software_versions.csv")
        run_metadata = os.path.join(self.output_dirs["report_directory"], "run_metadata.csv")

        freebayes_ncovtools_directory = os.path.join(self.output_dirs["report_directory"], "ncov_tools_freebayes")
        freebayes_metadata = os.path.join(freebayes_ncovtools_directory, "metadata.tsv")
        freebayes_ncovtools_data_directory = os.path.join(freebayes_ncovtools_directory, "data")
        freebayes_ncovtools_config = os.path.join(freebayes_ncovtools_directory, "config.yaml")

        job = concat_jobs(
            [
                bash.mkdir(freebayes_ncovtools_data_directory),
                Job(
                    input_files=[],
                    output_files=[
                        freebayes_readset_file_report,
                        freebayes_metadata
                    ],
                    command=f"""\\
head -n 1 {readset_file} > {freebayes_readset_file_report} && \\
echo -e "sample\\tct\\tdate" > {freebayes_metadata}"""
                )
            ]
        )

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            filtered_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")
            primer_trimmed_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.primerTrim.bam")
            consensus_directory = os.path.join(self.output_dirs["consensus_directory"], sample.name)
            freebayes_consensus = os.path.join(consensus_directory, f"{sample.name}.freebayes_calling.consensus.renamed.fasta")
            variant_directory = os.path.join(self.output_dirs["variants_directory"], sample.name)
            freebayes_variants = os.path.join(variant_directory, f"{sample.name}.freebayes_calling.consensus.vcf")

            freebayes_output_filtered_bam = os.path.join(freebayes_ncovtools_data_directory, os.path.basename(filtered_bam))
            freebayes_output_primer_trimmed_bam = os.path.join(freebayes_ncovtools_data_directory, os.path.basename(primer_trimmed_bam))
            freebayes_output_consensus = os.path.join(freebayes_ncovtools_data_directory, os.path.basename(freebayes_consensus))
            freebayes_output_variants = os.path.join(freebayes_ncovtools_data_directory, os.path.basename(freebayes_variants))

            job = concat_jobs(
                [
                    job,
                    Job(
                        input_files=[
                            filtered_bam,
                            primer_trimmed_bam,
                            freebayes_consensus,
                            freebayes_variants
                        ],
                        output_files=[
                            freebayes_output_filtered_bam,
                            freebayes_output_primer_trimmed_bam,
                            freebayes_output_consensus,
                            freebayes_output_variants
                        ],
                        command=f"""\
echo "Linking files for ncov_tools for sample {sample.name}..." && \\
if [ "$(ls -1 {filtered_bam})" != "" ] && [ "$(ls -1 {primer_trimmed_bam})" != "" ] && [ "$(ls -1 {freebayes_consensus})" != "" ] && [ "$(ls -1 {freebayes_variants})" != "" ];
  then
    ln -fs {os.path.relpath(filtered_bam, os.path.dirname(freebayes_output_filtered_bam))} {freebayes_output_filtered_bam} && \\
    ln -fs {os.path.relpath(primer_trimmed_bam, os.path.dirname(freebayes_output_primer_trimmed_bam))} {freebayes_output_primer_trimmed_bam} && \\
    ln -fs {os.path.relpath(freebayes_consensus, os.path.dirname(freebayes_output_consensus))} {freebayes_output_consensus} && \\
    ln -fs {os.path.relpath(freebayes_variants, os.path.dirname(freebayes_output_variants))} {freebayes_output_variants} && \\
    grep -w {sample.name} {readset_file} >> {freebayes_readset_file_report} && \\
    echo -e "{sample.name}\\tNA\\tNA" >> {freebayes_metadata}
fi"""
                    )
                ],
                samples=[sample]
            )

        jobs.append(
            concat_jobs(
                [
                    job,
                    ncovtools.run_ncovtools(
                        freebayes_output_filtered_bam,
                        freebayes_output_primer_trimmed_bam,
                        freebayes_output_consensus,
                        freebayes_output_variants,
                        readset_file,
                        freebayes_metadata,
                        freebayes_ncovtools_directory,
                        freebayes_ncovtools_config,
                        self.output_dir
                    ),
                    Job(
                        input_files=[],
                        output_files=[],
                        module_entries=[
                            ['prepare_report', 'module_R'],
                            ['prepare_report', 'module_CoVSeQ_tools'],
                            ['prepare_report', 'module_pandoc']
                        ],
                        command=f"""\\
module purge && \\
module load {global_conf.global_get('prepare_report', 'module_R')} {global_conf.global_get('prepare_report', 'module_CoVSeQ_tools')} {global_conf.global_get('prepare_report', 'module_pandoc')}"""
                    ),
                    covseq_tools.generate_report_tables(
                        freebayes_readset_file_report,
                        output_name_pattern=os.path.join("report", "report_metrics_freebayes")
                    ),
                    covseq_tools.render_report(
                        software_version,
                        run_metadata,
                        output_name_pattern=os.path.join("report", "report_metrics_freebayes"),
                        caller="freebayes"
                    )
                ],
                name="prepare_report." + global_conf.global_get('prepare_report', 'run_name', required=True)
            )
        )
        return jobs

    @property
    def step_list(self):
        return self.protocols()[self._protocol]

    def protocols(self):
        """
        Returns the protocol for the pipeline.
        Returns:
            dict: A dictionary of protocols for the pipeline.
        """
        return {'default': [
            self.host_reads_removal,
            self.kraken_analysis,
            self.cutadapt,
            self.mapping_bwa_mem_sambamba,
            self.sambamba_merge_sam_files,
            self.sambamba_filtering,
            self.ivar_trim_primers,
            self.covseq_metrics,
            self.freebayes_calling,
            self.ivar_calling,
            self.snpeff_annotate,
            self.ivar_create_consensus,
            self.bcftools_create_consensus,
            self.quast_consensus_metrics,
            self.rename_consensus_header_ivar,
            self.rename_consensus_header_freebayes,
            self.ncovtools_quickalign,
            self.prepare_table,
            self.prepare_report_ivar,
            self.prepare_report_freebayes,
            self.run_multiqc
        ]}

def main(parsed_args):
    """
    The function that will call this pipeline!
    Parameters:
        parsed_args (argparse.Namespace): The parsed arguments from the command line.
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

    pipeline = CoVSeq(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file, clean=clean, force=force, force_mem_per_cpu=force_mem_per_cpu, job_scheduler=job_scheduler, output_dir=output_dir, design_file=design_file, no_json=no_json, json_pt=json_pt, container=container)

    pipeline.submit_jobs()
