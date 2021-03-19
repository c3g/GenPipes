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

from bfx import bcftools
from bfx import bedtools
from bfx import bwa
from bfx import cutadapt
from bfx import fgbio
from bfx import freebayes
from bfx import gatk4
from bfx import htslib
from bfx import ivar
from bfx import kraken2
from bfx import multiqc
from bfx import quast
from bfx import sambamba
from bfx import samtools
from bfx import snpeff

from bfx import bash_cmd as bash

log = logging.getLogger(__name__)


class CoVSeQ(dnaseq.DnaSeqRaw):
    """
    CoVSeQ Pipeline
    ================

    pwet
    """

    def __init__(self, protocol=None):
        self._protocol = protocol
        # Add pipeline specific arguments
        super(CoVSeQ, self).__init__(protocol)


    def host_reads_removal(self):
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
            # trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            host_removal_directory = os.path.join("host_removal", readset.sample.name)
            readset_bam = os.path.join(host_removal_directory, readset.name + ".hybrid.sorted.bam")
            readset_bam_host_removed_sorted = os.path.join(host_removal_directory, readset.name + ".host_removed.sorted.bam")
            readset_bam_host_removed_sorted_index = os.path.join(host_removal_directory, readset.name + ".host_removed.sorted.bam.bai")
            readset_bam_host_removed_name_sorted = os.path.join(host_removal_directory, readset.name + ".host_removed.nsorted.bam")

            output_other = os.path.join(host_removal_directory, readset.name + ".host_removed.other.fastq.gz")
            output_single = os.path.join(host_removal_directory, readset.name + ".host_removed.single.fastq.gz")

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [
                    [readset.fastq1, readset.fastq2]
                    ]
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

                output_pair1 = os.path.join(host_removal_directory, readset.name + ".host_removed.pair1.fastq.gz")
                output_pair2 = os.path.join(host_removal_directory, readset.name + ".host_removed.pair2.fastq.gz")

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [
                    [readset.fastq1]
                    ]
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                output_pair1 = os.path.join(host_removal_directory, readset.name + ".host_removed.single.fastq.gz")
                output_pair2 = None

            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            jobs.append(
                concat_jobs([
                    bash.mkdir(host_removal_directory),
                    pipe_jobs([
                        bwa.mem(
                            fastq1,
                            fastq2,
                            read_group="'@RG" + \
                                "\\tID:" + readset.name + \
                                "\\tSM:" + readset.sample.name + \
                                "\\tLB:" + (readset.library if readset.library else readset.sample.name) + \
                                ("\\tPU:run" + readset.run + "_" + readset.lane if readset.run and readset.lane else "") + \
                                ("\\tCN:" + config.param('host_reads_removal', 'sequencing_center') if config.param('host_reads_removal', 'sequencing_center', required=False) else "") + \
                                ("\\tPL:" + config.param('host_reads_removal', 'sequencing_technology') if config.param('host_reads_removal', 'sequencing_technology', required=False) else "Illumina") + \
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
                            tmp_dir=config.param('host_reads_removal', 'tmp_dir', required=True),
                            other_options=config.param('host_reads_removal', 'sambamba_sort_other_options', required=False)
                            )
                        ]),
                    sambamba.view(
                            readset_bam,
                            readset_bam_host_removed_sorted,
                            options=config.param('host_reads_removal', 'sambamba_view_other_options')
                            ),
                    sambamba.sort(
                            readset_bam_host_removed_sorted,
                            readset_bam_host_removed_name_sorted,
                            tmp_dir=config.param('host_reads_removal', 'tmp_dir', required=True),
                            other_options=config.param('host_reads_removal', 'sambamba_name_sort_other_options', required=False)
                            ),
                    sambamba.index(
                        readset_bam_host_removed_sorted,
                        readset_bam_host_removed_sorted_index,
                        other_options=config.param('host_reads_removal', 'sambamba_index_other_options', required=False)
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
                    name="host_reads_removal." + readset.name,
                    samples=[readset.sample]
                    )
                )

        return jobs


    def kraken_analysis(self):
        """
        kraken
        """

        # TODO: include kraken analysis and output in metrics
        jobs = []
        for readset in self.readsets:
            host_removal_directory = os.path.join("host_removal", readset.sample.name)
            host_removal_file_prefix = os.path.join(host_removal_directory, readset.name)
            kraken_directory = os.path.join("metrics", "dna", readset.sample.name, "kraken_metrics")
            kraken_out_prefix = os.path.join(kraken_directory, readset.name)
            # readset_bam = os.path.join(host_removal_directory, readset.name, readset.name + ".nsorted.bam")
            # index_bam = os.path.join(host_removal_directory, readset.name, readset.name + ".nsorted.bam.bai")

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [
                    [host_removal_file_prefix + ".host_removed.pair1.fastq.gz", host_removal_file_prefix + ".host_removed.pair2.fastq.gz"]
                    ]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                unclassified_output = [kraken_out_prefix + ".unclassified_sequences_1.fastq", kraken_out_prefix + ".unclassified_sequences_2.fastq"]
                classified_output = [kraken_out_prefix + ".classified_sequences_1.fastq", kraken_out_prefix + ".classified_sequences_2.fastq"]

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [
                    [host_removal_file_prefix + ".host_removed.single.fastq.gz"]
                    ]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                unclassified_output = [kraken_out_prefix + ".unclassified_sequences.fastq"]
                classified_output = [kraken_out_prefix + ".classified_sequences.fastq"]

            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            jobs.append(
                concat_jobs([
                    bash.mkdir(kraken_directory),
                    kraken2.kraken2(
                        fastq1,
                        fastq2,
                        kraken_out_prefix,
                        other_options=config.param('kraken_analysis', 'kraken2_other_options'),
                        nthread=config.param('kraken_analysis', 'kraken2_threads'),
                        database=config.param('kraken_analysis', 'kraken2_database')
                        ),
                    Job(
                        input_files=unclassified_output + classified_output,
                        output_files=[s + ".gz" for s in unclassified_output + classified_output],
                        module_entries=[
                            ['pigz', 'module_pigz']
                        ],
                        command="""pigz -k -f -p {nthreads} {input_files}""".format(
                            input_files=" ".join(unclassified_output + classified_output),
                            nthreads=config.param('kraken_analysis', 'pigz_threads')
                            )
                        )
                    ],
                    name="kraken_analysis." + readset.name,
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
        """

        jobs = []

        for readset in self.readsets:
            host_removal_directory = os.path.join("host_removal", readset.sample.name)
            host_removal_file_prefix = os.path.join(host_removal_directory, readset.name)
            trim_directory = os.path.join("cleaned_raw_reads", readset.sample.name)
            trim_file_prefix = os.path.join(trim_directory, readset.name)

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [
                    [host_removal_file_prefix + ".host_removed.pair1.fastq.gz", host_removal_file_prefix + ".host_removed.pair2.fastq.gz"]
                    ]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                adapter_fwd = readset.adapter1
                adapter_rev = readset.adapter2

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [
                    [host_removal_file_prefix + ".host_removed.single.fastq.gz"]
                    ]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
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
            host_removal_directory = os.path.join("host_removal", readset.sample.name)
            host_removal_file_prefix = os.path.join(host_removal_directory, readset.name)
            trim_directory = os.path.join("cleaned_raw_reads", readset.sample.name)
            trim_file_prefix = os.path.join(trim_directory, readset.name)
            alignment_directory = os.path.join("alignment", readset.sample.name)
            readset_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam")
            index_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam.bai")

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [
                    [trim_file_prefix + ".trim.pair1.fastq.gz", trim_file_prefix + ".trim.pair2.fastq.gz"],
                    [host_removal_file_prefix + ".host_removed.pair1.fastq", host_removal_file_prefix + ".host_removed.pair2.fastq"]
                ]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [
                    [trim_file_prefix + ".trim.single.fastq.gz"],
                    [host_removal_file_prefix + ".host_removed.single.fastq"]
                ]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
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
                                ini_section='mapping_bwa_mem_sambamba'
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
            input_bams = [
                os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.bam"),
                os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam"),
                os.path.join(alignment_directory, sample.name + ".sorted.bam")
            ]

            for input_bam in input_bams:
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
                        name="sambamba_flagstat." + re.sub("\.bam$", "", os.path.basename(input_bam)),
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
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.bam")]
            ])

            output = os.path.join(alignment_directory, re.sub("\.bam$", ".BedGraph", os.path.basename(input_bam)))

            jobs.append(
                concat_jobs([
                    bash.mkdir(alignment_directory),
                    bedtools.genomecov(
                        input_bam,
                        output
                        )
                    ],
                    name="bedtools_genomecov." + sample.name,
                    samples=[sample]
                    )
                )

        return jobs


    def multiple_metrics_picard(self):
        # Calculates picard metrics from dnaseq pipeline but on raw AND on filtered bam file

        ##check the library status
        library = {}
        for readset in self.readsets:
            if not library.has_key(readset.sample):
                library[readset.sample] = "SINGLE_END"
            if readset.run_type == "PAIRED_END":
                library[readset.sample] = "PAIRED_END"

        jobs = []
        for sample in self.samples:
            picard_directory = os.path.join("metrics", "dna", sample.name, "picard_metrics")
            alignment_directory = os.path.join("alignment", sample.name)
            readset = sample.readsets[0]

            input_bams = [
                os.path.join(alignment_directory, sample.name + ".sorted.bam"),
                os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")
                ]
            mkdir_job = bash.mkdir(picard_directory, remove=True)

            for input in input_bams:
                jobs.append(
                    concat_jobs([
                        mkdir_job,
                        gatk4.collect_multiple_metrics(
                            input,
                            os.path.join(picard_directory, re.sub("\.bam$", "", os.path.basename(input)) + ".all.metrics"),
                            library_type=library[sample]
                            )
                        ],
                        name="multiple_metrics_raw_picard." + re.sub("\.bam$", "", os.path.basename(input)),
                        samples=[sample]
                        )
                    )

                jobs.append(
                    concat_jobs([
                        mkdir_job,
                        gatk4.collect_oxog_metrics(
                            input,
                            os.path.join(picard_directory, re.sub("\.bam$", "", os.path.basename(input)) + ".oxog_metrics.txt")
                            )
                        ],
                        name="picard_collect_oxog_metrics." + re.sub("\.bam$", "", os.path.basename(input)),
                        samples=[sample]
                        )
                    )

                jobs.append(
                    concat_jobs([
                        mkdir_job,
                        gatk4.collect_gcbias_metrics(
                            input,
                            os.path.join(picard_directory, re.sub("\.bam$", "", os.path.basename(input)) + ".qcbias_metrics.txt"),
                            os.path.join(picard_directory, re.sub("\.bam$", "", os.path.basename(input)) + ".qcbias_metrics.pdf"),
                            os.path.join(picard_directory, re.sub("\.bam$", "", os.path.basename(input)) + ".qcbias_summary_metrics.txt")
                            )
                        ],
                        name="picard_collect_gcbias_metrics." + re.sub("\.bam$", "", os.path.basename(input)),
                        samples=[sample]
                        )
                    )

        return jobs



    def covseq_metrics(self):
        """
        Multiple metrcis from dnaseq
        """

        jobs = []

        jobs.extend(self.multiple_metrics_picard())
        # jobs.extend(self.metrics_dna_picard_metrics())
        jobs.extend(self.metrics_dna_sample_qualimap())
        jobs.extend(self.metrics_sambamba_flagstat())
        jobs.extend(self.metrics_bedtools_genomecov())
        # jobs.extend(self.metrics_dna_sambamba_flagstat())
        jobs.extend(self.picard_calculate_hs_metrics())
        jobs.extend(self.metrics())

        return jobs


    def ivar_calling(self):
        """
        ivar calling
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            variant_directory = os.path.join("variant", sample.name)

            [input_bam] = self.select_input_files([
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.bam")]
            ])

            output_prefix = os.path.join(variant_directory, sample.name) + ".variants"
            output_tsv = output_prefix + ".tsv"
            output_vcf = os.path.join(variant_directory, re.sub("\.bam$", "", os.path.basename(input_bam))) + ".vcf"

            jobs.append(
                concat_jobs([
                    bash.mkdir(variant_directory),
                    pipe_jobs([
                        samtools.mpileup(
                            input_bam,
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
                    name="ivar_call_variants." + sample.name,
                    samples=[sample],
                    removable_files=[output_vcf]
                    )
                )

        return jobs


    def freebayes_calling(self):
        """
        freebayes calling
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            variant_directory = os.path.join("variant", sample.name)

            [input_bam] = self.select_input_files([
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.bam")]
            ])

            output_prefix = os.path.join(variant_directory, sample.name) + ".freebayes_calling"
            output_gvcf = output_prefix + ".gvcf"

            output_masks = output_prefix + ".mask.txt"
            output_variants = output_prefix + ".variants.vcf"
            output_variants_norm = re.sub(r"\.vcf$", ".norm.vcf", output_variants)
            output_consensus = output_prefix + ".consensus.vcf"
            output_consensus_norm = re.sub(r"\.vcf$", ".norm.vcf", output_consensus)
            output_ambiguous_norm = output_prefix + ".ambiguous.norm.vcf"
            output_fixed_norm = output_prefix + ".fixed.norm.vcf"

            jobs.append(
                concat_jobs([
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
                        ),
                    Job(
                        input_files=[output_variants, output_consensus],
                        output_files=[output_variants_norm, output_variants_norm + ".gz", output_consensus_norm, output_consensus_norm + ".gz", output_ambiguous_norm, output_ambiguous_norm + ".gz", output_fixed_norm, output_fixed_norm + ".gz"],
                        module_entries=[
                            [freebayes, 'module_bcftools'],
                            [freebayes, 'module_htslib']
                        ],
                        command="""\\
for v in "variants" "consensus"
do
    bcftools norm -f {reference_genome} {output_prefix}.$v.vcf > {output_prefix}.$v.norm.vcf
    bgzip -f {output_prefix}.$v.norm.vcf
    tabix -f -p vcf {output_prefix}.$v.norm.vcf.gz
done && \\
for vt in "ambiguous" "fixed"
do
    cat {output_prefix}.consensus.norm.vcf | awk -v vartag=ConsensusTag=$vt '$0 ~ /^#/ || $0 ~ vartag' > {output_prefix}.$vt.norm.vcf
    bgzip -f {output_prefix}.$vt.norm.vcf
    tabix -f -p vcf {output_prefix}.$vt.norm.vcf.gz
done""".format(
    reference_genome=config.param("DEFAULT", 'genome_fasta', type='filepath'),
    output_prefix=output_prefix
    )
            )
                    ],
                    name="freebayes_call_variants." + sample.name,
                    samples=[sample],
                    removable_files=[output_gvcf]
                    )
                )

        return jobs


    def snpeff_annotate(self):
        """
        Consensus annotation with SnpEff
        """

        jobs = []

        for sample in self.samples:
            variant_directory = os.path.join("variant", sample.name)
            metrics_prefix = os.path.join("metrics", "dna", sample.name, "snpeff_metrics", sample.name + ".snpEff")

            [input_vcf] = self.select_input_files([
                [os.path.join(variant_directory, sample.name + ".sorted.filtered.primerTrim.vcf.gz")],
                [os.path.join(variant_directory, sample.name + ".sorted.filtered.vcf.gz")],
                [os.path.join(variant_directory, sample.name + ".sorted.vcf.gz")]
            ])

            output_vcf = os.path.join(variant_directory, re.sub("\.vcf.gz$", ".annotate.vcf", os.path.basename(input_vcf)))

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(metrics_prefix)),
                    snpeff.snpeff_annotate(
                        input_vcf,
                        output_vcf,
                        metrics_prefix
                        ),
                    htslib.bgzip_tabix(
                        output_vcf,
                        output_vcf + ".gz"
                        )
                    ],
                    name="snpeff_annotate." + sample.name,
                    samples=[sample],
                    removable_files=[output_vcf]
                    )
                )

        return jobs

    def ivar_create_consensus(self):
        """
        Create consensus with ivar through a samtools mpileup
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            consensus_directory = os.path.join("consensus", sample.name)
            [input_bam] = self.select_input_files([
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.bam")]
            ])

            output_prefix = os.path.join(consensus_directory, re.sub("\.bam$", ".consensus", os.path.basename(input_bam)))

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(output_prefix)),
                    pipe_jobs([
                        samtools.mpileup(
                            input_bam,
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
                    name="ivar_create_consensus." + sample.name,
                    samples=[sample],
                    removable_files=[output_prefix + ".fa"]
                    )
                )

        return jobs

    def bcftools_create_consensus(self):
        """
        bcftools consensus creation
        """

        jobs = []
        for sample in self.samples:
            consensus_directory = os.path.join("consensus", sample.name)
            variant_directory = os.path.join("variant", sample.name)


            input_prefix = os.path.join(variant_directory, sample.name) + ".freebayes_calling"

            intput_masks = input_prefix + ".mask.txt"
            input_ambiguous_norm = input_prefix + ".ambiguous.norm.vcf.gz"
            input_fixed_norm = input_prefix + ".fixed.norm.vcf.gz"

            output_prefix = os.path.join(consensus_directory, sample.name) + ".freebayes_calling"
            output_ambiguous_fasta = output_prefix + ".ambiguous.fasta"
            output_consensus_fasta = output_prefix + ".consensus.fasta"

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(output_prefix)),
                    bcftools.consensus(
                        input_ambiguous_norm,
                        output_ambiguous_fasta,
                        "-f " + config.param("DEFAULT", 'genome_fasta', type='filepath') + " -I "
                        ),
                    pipe_jobs([
                        bcftools.consensus(
                            input_fixed_norm,
                            None,
                            "-f " + output_ambiguous_fasta + " -m " + intput_masks
                            ),
                        Job(
                            input_files=[],
                            output_files=[output_consensus_fasta],
                            module_entries=[],
                            command="""\\
sed s/{reference_genome_name}/{sample_name}/ > {output_consensus_fasta}""".format(
    reference_genome_name=config.param("DEFAULT", 'assembly_synonyms'),
    sample_name=sample.name,
    output_consensus_fasta=output_consensus_fasta
    )
                            )
                        ])
                    ],
                    name="bcftools_create_consensus." + sample.name,
                    samples=[sample]
                    )
                )

        return jobs


    def quast_consensus_metrics(self):
        """
        Generate QUAST metrics on consensus
        """

        jobs = []
        for sample in self.samples:
            consensus_directory = os.path.join("consensus", sample.name)
            output_dir = os.path.join("metrics", "dna", sample.name, "quast_metrics")
            [input_fa] = self.select_input_files([
                [os.path.join(consensus_directory, sample.name + ".sorted.filtered.primerTrim.consensus.fa")],
                [os.path.join(consensus_directory, sample.name + ".sorted.filtered.consensus.fa")],
                [os.path.join(consensus_directory, sample.name + ".sorted.consensus.fa")]
            ])

            jobs.append(
                concat_jobs([
                    bash.mkdir(output_dir),
                    quast.quast(
                        input_fa,
                        output_dir,
                        prefix="quast_consensus_metrics"
                        )
                    ],
                    name="quast_consensus_metrics." + sample.name,
                    samples=[sample]
                    )
                )

        return jobs


    def rename_consensus_header(self):
        """
        Rename reads headers
        """

        jobs = []

        job = bash.mkdir(os.path.join("consensus"))
        for sample in self.samples:
            consensus_directory = os.path.join("consensus", sample.name)
            [input_fa] = self.select_input_files([
                [os.path.join(consensus_directory, sample.name + ".sorted.filtered.primerTrim.consensus.fa")],
                [os.path.join(consensus_directory, sample.name + ".sorted.filtered.consensus.fa")],
                [os.path.join(consensus_directory, sample.name + ".sorted.consensus.fa")]
            ])
            quast_directory = os.path.join("metrics", "dna", sample.name, "quast_metrics")
            quast_html = os.path.join(quast_directory, "report.html")
            quast_tsv = os.path.join(quast_directory, "report.tsv")

            output_fa = os.path.join(consensus_directory, sample.name + ".consensus.fasta")
            output_status_fa = os.path.join(consensus_directory, """{sample_name}.consensus.{technology}.{status}.fasta""".format(sample_name=sample.name, technology=config.param('rename_consensus_header', 'seq_technology', required=False), status="${STATUS}"))

            variant_directory = os.path.join("variant", sample.name)
            [input_vcf] = self.select_input_files([
                [os.path.join(variant_directory, sample.name + ".sorted.filtered.primerTrim.vcf.gz")],
                [os.path.join(variant_directory, sample.name + ".sorted.filtered.vcf.gz")],
                [os.path.join(variant_directory, sample.name + ".sorted.vcf.gz")]
            ])
            annotated_vcf = os.path.join(variant_directory, re.sub("\.vcf.gz$", ".annotate.vcf", os.path.basename(input_vcf)))

            alignment_directory = os.path.join("alignment", sample.name)
            [input_bam] = self.select_input_files([
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.bam")]
            ])

            bedgraph_file = os.path.join(alignment_directory, re.sub("\.bam$", ".BedGraph", os.path.basename(input_bam)))

            job = concat_jobs([
                job,
                bash.mkdir(os.path.dirname(output_fa)),
                Job(
                    input_files=[quast_tsv, quast_html],
                    output_files=[],
                    command="""\\
cons_len=`grep -oP "Total length \(>= 0 bp\)\\t\K.*?(?=$)" {quast_tsv}`
N_count=`grep -oP "# N's\\",\\"quality\\":\\"Less is better\\",\\"values\\":\[\K.*?(?=])" {quast_html}`
cons_perc_N=`echo "scale=2; 100*$N_count/$cons_len" | bc -l`
frameshift=`if grep -q "frameshift_variant" {annotated_vcf}; then echo "FLAG"; fi`
genome_size=`awk '{{print $2}}' {genome_file}`
bam_cov50X=`awk '{{if ($4 > 50) {{count = count + $3-$2}}}} END {{if (count) {{print count}} else {{print 0}}}}' {bedgraph_file}`
bam_cov50X=`echo "scale=2; 100*$bam_cov50X/$genome_size" | bc -l`
STATUS=`awk -v bam_cov50X=$bam_cov50X -v frameshift=$frameshift -v cons_perc_N=$cons_perc_N 'BEGIN {{ if (cons_perc_N < 1 && frameshift != "FLAG" && bam_cov50X >= 90) {{print "pass"}}  else if (cons_perc_N > 5) {{print "rej"}} else if ((cons_perc_N >= 1 && cons_perc_N <= 5) || frameshift == "FLAG" || bam_cov50X < 90) {{print "flag"}} }}'`
export STATUS""".format(
    quast_html=quast_html,
    quast_tsv=quast_tsv,
    genome_file=config.param('DEFAULT', 'igv_genome', required=False),
    annotated_vcf=annotated_vcf,
    bedgraph_file=bedgraph_file
    )
                    ),
                Job(
                    input_files=[input_fa],
                    output_files=[output_fa],
                    command="""\\
awk '/^>/{{print ">{country}/{province}-{sample}/{year} seq_method:{seq_method}|assemb_method:{assemb_method}|snv_call_method:{snv_call_method}"; next}}{{print}}' < {input_fa} > {output_status_fa} && \\
ln -sf {output_status_fa_basename} {output_fa}
""".format(
    country=config.param('rename_consensus_header', 'country', required=False),
    province=config.param('rename_consensus_header', 'province', required=False),
    year=config.param('rename_consensus_header', 'year', required=False),
    seq_method=config.param('rename_consensus_header', 'seq_method', required=False),
    assemb_method=config.param('rename_consensus_header', 'assemb_method', required=False),
    snv_call_method=config.param('rename_consensus_header', 'snv_call_method', required=False),
    sample=sample.name,
    input_fa=input_fa,
    output_status_fa_basename=os.path.basename(output_status_fa),
    output_status_fa=output_status_fa,
    output_fa=output_fa
    )
                    )
            ],
            name="rename_consensus_header",
            samples=[sample]
            )

        jobs.append(job)

#             jobs.append(
#                 concat_jobs([
#                     bash.mkdir(os.path.dirname(output_fa)),
#                     Job(
#                         input_files=[quast_tsv, quast_html],
#                         output_files=[],
#                         command="""\\
# cons_len=`grep -oP "Total length \(>= 0 bp\)\\t\K.*?(?=$)" {quast_tsv}`
# N_count=`grep -oP "# N's\\",\\"quality\\":\\"Less is better\\",\\"values\\":\[\K.*?(?=])" {quast_html}`
# cons_perc_N=`echo "scale=2; 100*$N_count/$cons_len" | bc -l`
# frameshift=`if grep -q "frameshift_variant" {annotated_vcf}; then echo "FLAG"; fi`
# genome_size=`awk '{{print $2}}' {genome_file}`
# bam_cov50X=`awk '{{if ($4 > 50) {{count = count + $3-$2}}}} END {{if (count) {{print count}} else {{print 0}}}}' {bedgraph_file}`
# bam_cov50X=`echo "scale=2; 100*$bam_cov50X/$genome_size" | bc -l`
# STATUS=`awk -v bam_cov50X=$bam_cov50X -v frameshift=$frameshift -v cons_perc_N=$cons_perc_N 'BEGIN {{ if (cons_perc_N < 1 && frameshift != "FLAG" && bam_cov50X >= 90) {{print "pass"}}  else if (cons_perc_N > 5) {{print "rej"}} else if ((cons_perc_N >= 1 && cons_perc_N <= 5) || frameshift == "FLAG" || bam_cov50X < 90) {{print "flag"}} }}'`
# export STATUS""".format(
#     quast_html=quast_html,
#     quast_tsv=quast_tsv,
#     genome_file=config.param('DEFAULT', 'igv_genome', required=False),
#     annotated_vcf=annotated_vcf,
#     bedgraph_file=bedgraph_file
#     )
#                         ),
#                     Job(
#                         input_files=[input_fa],
#                         output_files=[output_fa],
#                         command="""\\
# awk '/^>/{{print ">{country}/{province}-{sample}/{year} seq_method:{seq_method}|assemb_method:{assemb_method}|snv_call_method:{snv_call_method}"; next}}{{print}}' < {input_fa} > {output_status_fa} && \\
# ln -sf {output_status_fa_basename} {output_fa}
# """.format(
#     country=config.param('rename_consensus_header', 'country', required=False),
#     province=config.param('rename_consensus_header', 'province', required=False),
#     year=config.param('rename_consensus_header', 'year', required=False),
#     seq_method=config.param('rename_consensus_header', 'seq_method', required=False),
#     assemb_method=config.param('rename_consensus_header', 'assemb_method', required=False),
#     snv_call_method=config.param('rename_consensus_header', 'snv_call_method', required=False),
#     sample=sample.name,
#     input_fa=input_fa,
#     output_status_fa_basename=os.path.basename(output_status_fa),
#     output_status_fa=output_status_fa,
#     output_fa=output_fa
#     )
#                         )
#                 ],
#                 name="rename_consensus_header." + sample.name,
#                 samples=[sample]
#                 )
#             )

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
            self.host_reads_removal,
            self.kraken_analysis,
            self.cutadapt,
            self.mapping_bwa_mem_sambamba,
            self.sambamba_merge_sam_files,
            self.sambamba_filtering,
            # self.fgbio_trim_primers,
            self.ivar_trim_primers,
            self.covseq_metrics,
            self.freebayes_calling,
            self.ivar_calling,
            self.snpeff_annotate,
            self.ivar_create_consensus,
            self.bcftools_create_consensus,
            self.quast_consensus_metrics,
            self.rename_consensus_header
            # self.run_multiqc
        ]

if __name__ == '__main__':
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        CoVSeQ()
