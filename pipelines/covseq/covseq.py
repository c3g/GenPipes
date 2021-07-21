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
import subprocess

from pipelines.dnaseq import dnaseq

from bfx import bcftools
from bfx import bedtools
from bfx import bwa
from bfx import covseq_tools
from bfx import cutadapt
from bfx import fgbio
from bfx import freebayes
from bfx import gatk4
from bfx import htslib
from bfx import ivar
from bfx import kraken2
from bfx import multiqc
from bfx import ncovtools
from bfx import quast
from bfx import sambamba
from bfx import samtools
from bfx import snpeff

from bfx import bash_cmd as bash
from core.config import config

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
            # output_variants_norm = re.sub(r"\.vcf$", ".norm.vcf.gz", output_variants)
            output_consensus = output_prefix + ".consensus.vcf"
            output_consensus_norm = re.sub(r"\.vcf$", ".norm.vcf.gz", output_consensus)
            # output_ambiguous_norm = output_prefix + ".ambiguous.norm.vcf.gz"
            # output_fixed_norm = output_prefix + ".fixed.norm.vcf.gz"


            job = concat_jobs([
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
                ])

            for file in [output_variants, output_consensus]:
                output_vcf_gz = re.sub(r"\.vcf$", ".norm.vcf.gz", file)
                job = concat_jobs([
                    job,
                    pipe_jobs([
                        bcftools.norm(
                            file,
                            None,
                            "-f " + config.param("DEFAULT", 'genome_fasta', type='filepath'),
                            ini_section='freebayes_call_variants'
                            ),
                        htslib.bgzip_tabix(
                            None,
                            output_vcf_gz
                            )
                        ])
                    ])

            for flag in ["ambiguous", "fixed"]:
                output_vcf_gz = re.sub(r"\.consensus\.norm\.vcf\.gz$", ".%s.norm.vcf.gz" % (flag), output_consensus_norm)
                job = concat_jobs([
                    job,
                    pipe_jobs([
                        bash.cat(
                            output_consensus_norm,
                            None,
                            zip=True
                        ),
                        bash.awk(
                            None,
                            None,
                            "-v vartag=ConsensusTag=%s '$0 ~ /^#/ || $0 ~ vartag'" % (flag)),
                        htslib.bgzip_tabix(
                            None,
                            output_vcf_gz
                            )
                        ])
                    ])

            job.name = "freebayes_call_variants." + sample.name
            job.samples = [sample]
            job.removable_files = [output_gvcf]
            jobs.append(job)

        return jobs


    def snpeff_annotate(self):
        """
        Consensus annotation with SnpEff
        """

        jobs = []

        for sample in self.samples:
            variant_directory = os.path.join("variant", sample.name)
            metrics_ivar_prefix = os.path.join("metrics", "dna", sample.name, "snpeff_metrics_ivar", sample.name + ".snpEff")
            metrics_freebayes_prefix = os.path.join("metrics", "dna", sample.name, "snpeff_metrics_freebayes", sample.name + ".snpEff")

            [ivar_vcf] = self.select_input_files([
                [os.path.join(variant_directory, sample.name + ".sorted.filtered.primerTrim.vcf.gz")],
                [os.path.join(variant_directory, sample.name + ".sorted.filtered.vcf.gz")],
                [os.path.join(variant_directory, sample.name + ".sorted.vcf.gz")]
            ])

            ivar_output_vcf = os.path.join(variant_directory, re.sub("\.vcf.gz$", ".annotate.vcf", os.path.basename(ivar_vcf)))

            freebayes_vcf = os.path.join(variant_directory, sample.name) + ".freebayes_calling.fixed.norm.vcf.gz"
            freebayes_output_vcf = os.path.join(variant_directory, re.sub("\.vcf.gz$", ".annotate.vcf", os.path.basename(freebayes_vcf)))

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(metrics_ivar_prefix)),
                    snpeff.snpeff_annotate(
                        ivar_vcf,
                        ivar_output_vcf,
                        metrics_ivar_prefix
                        ),
                    htslib.bgzip_tabix(
                        ivar_output_vcf,
                        ivar_output_vcf + ".gz"
                        ),
                    bash.mkdir(os.path.dirname(metrics_freebayes_prefix)),
                    snpeff.snpeff_annotate(
                        freebayes_vcf,
                        freebayes_output_vcf,
                        metrics_freebayes_prefix
                        ),
                    htslib.bgzip_tabix(
                        freebayes_output_vcf,
                        freebayes_output_vcf + ".gz"
                        )
                    ],
                    name="snpeff_annotate." + sample.name,
                    samples=[sample]
                    # removable_files=[output_vcf]
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
            ivar_output_dir = os.path.join("metrics", "dna", sample.name, "quast_metrics_ivar")
            [ivar_consensus] = self.select_input_files([
                [os.path.join(consensus_directory, sample.name + ".sorted.filtered.primerTrim.consensus.fa")],
                [os.path.join(consensus_directory, sample.name + ".sorted.filtered.consensus.fa")],
                [os.path.join(consensus_directory, sample.name + ".sorted.consensus.fa")]
            ])
            freebayes_output_dir = os.path.join("metrics", "dna", sample.name, "quast_metrics_freebayes")
            freebayes_consensus = os.path.join(consensus_directory, sample.name) + ".freebayes_calling.consensus.fasta"

            jobs.append(
                concat_jobs([
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
                    name="quast_consensus_metrics." + sample.name,
                    samples=[sample]
                    )
                )

        return jobs


    def rename_consensus_header_ivar(self):
        """
        Rename reads headers
        """

        jobs = []

        for sample in self.samples:
            consensus_directory = os.path.join("consensus", sample.name)
            [ivar_consensus] = self.select_input_files([
                [os.path.join(consensus_directory, sample.name + ".sorted.filtered.primerTrim.consensus.fa")],
                [os.path.join(consensus_directory, sample.name + ".sorted.filtered.consensus.fa")],
                [os.path.join(consensus_directory, sample.name + ".sorted.consensus.fa")]
            ])
            quast_ivar_directory = os.path.join("metrics", "dna", sample.name, "quast_metrics_ivar")
            quast_ivar_html = os.path.join(quast_ivar_directory, "report.html")
            quast_ivar_tsv = os.path.join(quast_ivar_directory, "report.tsv")

            ivar_output_fa = os.path.join(consensus_directory, sample.name + ".consensus.fasta")
            ivar_output_status_fa = os.path.join(consensus_directory, """{sample_name}.consensus.{technology}.{status}.fasta""".format(sample_name=sample.name, technology=config.param('rename_consensus_header', 'sequencing_technology', required=False), status="${IVAR_STATUS}"))

            variant_directory = os.path.join("variant", sample.name)
            [ivar_vcf] = self.select_input_files([
                [os.path.join(variant_directory, sample.name + ".sorted.filtered.primerTrim.vcf.gz")],
                [os.path.join(variant_directory, sample.name + ".sorted.filtered.vcf.gz")],
                [os.path.join(variant_directory, sample.name + ".sorted.vcf.gz")]
            ])
            ivar_annotated_vcf = os.path.join(variant_directory, re.sub("\.vcf.gz$", ".annotate.vcf", os.path.basename(ivar_vcf)))

            alignment_directory = os.path.join("alignment", sample.name)
            [input_bam] = self.select_input_files([
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.bam")]
            ])

            bedgraph_file = os.path.join(alignment_directory, re.sub("\.bam$", ".BedGraph", os.path.basename(input_bam)))

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(ivar_output_fa)),
                    Job(
                        input_files=[quast_ivar_tsv, quast_ivar_html],
                        output_files=[],
                        command="""\\
ivar_cons_len=`grep -oP "Total length \(>= 0 bp\)\\t\K.*?(?=$)" {quast_ivar_tsv}`
ivar_N_count=`grep -oP "# N's\\",\\"quality\\":\\"Less is better\\",\\"values\\":\[\K.*?(?=])" {quast_ivar_html}`
ivar_cons_perc_N=`echo "scale=2; 100*$ivar_N_count/$ivar_cons_len" | bc -l`
ivar_frameshift=`if grep -q "frameshift_variant" {ivar_annotated_vcf}; then echo "FLAG"; fi`
genome_size=`awk '{{print $2}}' {genome_file}`
bam_cov50X=`awk '{{if ($4 > 50) {{count = count + $3-$2}}}} END {{if (count) {{print count}} else {{print 0}}}}' {bedgraph_file}`
bam_cov50X=`echo "scale=2; 100*$bam_cov50X/$genome_size" | bc -l`
IVAR_STATUS=`awk -v bam_cov50X=$bam_cov50X -v ivar_frameshift=$ivar_frameshift -v ivar_cons_perc_N=$ivar_cons_perc_N 'BEGIN {{ if (ivar_cons_perc_N < 1 && ivar_frameshift != "FLAG" && bam_cov50X >= 90) {{print "pass"}}  else if (ivar_cons_perc_N > 5) {{print "rej"}} else if ((ivar_cons_perc_N >= 1 && ivar_cons_perc_N <= 5) || ivar_frameshift == "FLAG" || bam_cov50X < 90) {{print "flag"}} }}'`
export IVAR_STATUS""".format(
    quast_ivar_html=quast_ivar_html,
    quast_ivar_tsv=quast_ivar_tsv,
    genome_file=config.param('DEFAULT', 'igv_genome', required=False),
    ivar_annotated_vcf=ivar_annotated_vcf,
    bedgraph_file=bedgraph_file
    )
                        ),
                    Job(
                        input_files=[ivar_consensus],
                        output_files=[ivar_output_fa],
                        command="""\\
awk '/^>/{{print ">{country}/{province}-{sample}/{year} seq_method:{seq_method}|assemb_method:ivar|snv_call_method:ivar"; next}}{{print}}' < {ivar_consensus} > {ivar_output_status_fa} && \\
ln -sf {ivar_output_status_fa_basename} {ivar_output_fa}""".format(
    country=config.param('rename_consensus_header', 'country', required=False),
    province=config.param('rename_consensus_header', 'province', required=False),
    year=config.param('rename_consensus_header', 'year', required=False),
    seq_method=config.param('rename_consensus_header', 'seq_method', required=False),
    sample=sample.name,
    ivar_consensus=ivar_consensus,
    ivar_output_status_fa_basename=os.path.basename(ivar_output_status_fa),
    ivar_output_status_fa=ivar_output_status_fa,
    ivar_output_fa=ivar_output_fa
    )
                        )
                ],
                name="rename_consensus_header." + sample.name,
                samples=[sample]
                )
            )

        return jobs


    def rename_consensus_header_freebayes(self):
        """
        Rename reads headers
        """

        jobs = []

        for sample in self.samples:
            consensus_directory = os.path.join("consensus", sample.name)
            freebayes_consensus = os.path.join(consensus_directory, sample.name) + ".freebayes_calling.consensus.fasta"
            
            quast_freebayes_directory = os.path.join("metrics", "dna", sample.name, "quast_metrics_freebayes")
            quast_freebayes_html = os.path.join(quast_freebayes_directory, "report.html")
            quast_freebayes_tsv = os.path.join(quast_freebayes_directory, "report.tsv")

            freebayes_output_fa = os.path.join(consensus_directory, sample.name) + ".freebayes_calling.consensus.renamed.fasta"
            freebayes_output_status_fa = os.path.join(consensus_directory, """{sample_name}.freebayes_calling.consensus.{technology}.{status}.fasta""".format(sample_name=sample.name, technology=config.param('rename_consensus_header', 'sequencing_technology', required=False), status="${FREEBAYES_STATUS}"))

            variant_directory = os.path.join("variant", sample.name)
            freebayes_vcf = os.path.join(variant_directory, sample.name) + ".freebayes_calling.fixed.norm.vcf.gz"
            freebayes_annotated_vcf = os.path.join(variant_directory, re.sub("\.vcf.gz$", ".annotate.vcf", os.path.basename(freebayes_vcf)))

            alignment_directory = os.path.join("alignment", sample.name)
            [input_bam] = self.select_input_files([
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.bam")]
            ])

            bedgraph_file = os.path.join(alignment_directory, re.sub("\.bam$", ".BedGraph", os.path.basename(input_bam)))

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(freebayes_output_fa)),
                    Job(
                        input_files=[quast_freebayes_tsv, quast_freebayes_html],
                        output_files=[],
                        command="""\\
freebayes_cons_len=`grep -oP "Total length \(>= 0 bp\)\\t\K.*?(?=$)" {quast_freebayes_tsv}`
freebayes_N_count=`grep -oP "# N's\\",\\"quality\\":\\"Less is better\\",\\"values\\":\[\K.*?(?=])" {quast_freebayes_html}`
freebayes_cons_perc_N=`echo "scale=2; 100*$freebayes_N_count/$freebayes_cons_len" | bc -l`
freebayes_frameshift=`if grep -q "frameshift_variant" {freebayes_annotated_vcf}; then echo "FLAG"; fi`
genome_size=`awk '{{print $2}}' {genome_file}`
bam_cov50X=`awk '{{if ($4 > 50) {{count = count + $3-$2}}}} END {{if (count) {{print count}} else {{print 0}}}}' {bedgraph_file}`
bam_cov50X=`echo "scale=2; 100*$bam_cov50X/$genome_size" | bc -l`
FREEBAYES_STATUS=`awk -v bam_cov50X=$bam_cov50X -v freebayes_frameshift=$freebayes_frameshift -v freebayes_cons_perc_N=$freebayes_cons_perc_N 'BEGIN {{ if (freebayes_cons_perc_N < 1 && freebayes_frameshift != "FLAG" && bam_cov50X >= 90) {{print "pass"}}  else if (freebayes_cons_perc_N > 5) {{print "rej"}} else if ((freebayes_cons_perc_N >= 1 && freebayes_cons_perc_N <= 5) || freebayes_frameshift == "FLAG" || bam_cov50X < 90) {{print "flag"}} }}'`
export FREEBAYES_STATUS""".format(
    quast_freebayes_html=quast_freebayes_html,
    quast_freebayes_tsv=quast_freebayes_tsv,
    genome_file=config.param('DEFAULT', 'igv_genome', required=False),
    freebayes_annotated_vcf=freebayes_annotated_vcf,
    bedgraph_file=bedgraph_file
    )
                        ),
                    Job(
                        input_files=[freebayes_consensus],
                        output_files=[freebayes_output_fa],
                        command="""\\
awk '/^>/{{print ">{country}/{province}-{sample}/{year} seq_method:{seq_method}|assemb_method:bcftools|snv_call_method:freebayes"; next}}{{print}}' < {freebayes_consensus} > {freebayes_output_status_fa} && \\
ln -sf {freebayes_output_status_fa_basename} {freebayes_output_fa}""".format(
    country=config.param('rename_consensus_header', 'country', required=False),
    province=config.param('rename_consensus_header', 'province', required=False),
    year=config.param('rename_consensus_header', 'year', required=False),
    seq_method=config.param('rename_consensus_header', 'seq_method', required=False),
    sample=sample.name,
    freebayes_consensus=freebayes_consensus,
    freebayes_output_status_fa_basename=os.path.basename(freebayes_output_status_fa),
    freebayes_output_status_fa=freebayes_output_status_fa,
    freebayes_output_fa=freebayes_output_fa
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

            input_dep += [
                input_oxog,
                input_qcbias,
                input_all_picard
            ]

            inputs += [os.path.join(metrics_directory, sample.name)]

        output = os.path.join(metrics_directory, "multiqc_report")

        job = multiqc.run(
            input_dep,
            output
            )
        job.name = "multiqc_all_samples"
        job.samples = self.samples

        jobs.append(job)

        return jobs


    def ncovtools_quickalign(self):

        jobs = []

        for sample in self.samples:
            ncovtools_quickalign_directory = os.path.join("metrics", "dna", sample.name, "ncovtools_quickalign")
            output = os.path.join(ncovtools_quickalign_directory, sample.name + "_ivar_vs_freebayes.vcf")
            consensus_directory = os.path.join("consensus", sample.name)
            ivar_consensus = os.path.join(consensus_directory, sample.name + ".consensus.fasta")
            freebayes_consensus = os.path.join(consensus_directory, sample.name) + ".freebayes_calling.consensus.renamed.fasta"
            jobs.append(
                    concat_jobs([
                        bash.mkdir(ncovtools_quickalign_directory),
                        Job(
                            input_files=[ivar_consensus, freebayes_consensus],
                            output_files=[output],
                            module_entries=[
                                ['ncovtools_quickalign', 'module_python'],
                                ['ncovtools_quickalign', 'module_ncov_random_scripts'],
                            ],
                            command="""\\
quick_align.py -r {ivar_consensus} -g {freebayes_consensus} -o vcf > {output}""".format(
    ivar_consensus=ivar_consensus,
    freebayes_consensus=freebayes_consensus,
    output=output
    )
                            )
                    ],
                    name="ncovtools_quickalign." + sample.name,
                    samples=[sample]
                    )
                )

        return jobs

    def prepare_table(self):

        jobs = []

        metrics_directory = os.path.join("metrics", "dna")

        readset_file=os.path.relpath(self.args.readsets.name, self.output_dir)

        run_metadata = os.path.join("report", "run_metadata.csv")
        software_version = os.path.join("report", "software_versions.csv")

        modules = []
        # Retrieve all unique module version values in config files
        # assuming that all module key names start with "module_"
        for section in config.sections():
            for name, value in config.items(section):
                if re.search("^module_", name) and value not in modules:
                    modules.append(value)

        # Finding all kraken outputs
        kraken_outputs = []
        library = {}
        for readset in self.readsets:
            ##check the library status
            if not library.has_key(readset.sample):
                library[readset.sample] = "SINGLE_END"
            if readset.run_type == "PAIRED_END":
                library[readset.sample] = "PAIRED_END"

            kraken_directory = os.path.join("metrics", "dna", readset.sample.name, "kraken_metrics")
            kraken_out_prefix = os.path.join(kraken_directory, readset.name)
            kraken_output = kraken_out_prefix + ".kraken2_output"
            kraken_report = kraken_out_prefix + ".kraken2_report"
            kraken_outputs.extend((kraken_output, kraken_report))
            if readset.run_type == "PAIRED_END":
                unclassified_output_1 = kraken_out_prefix + ".unclassified_sequences_1.fastq"
                unclassified_output_2 = kraken_out_prefix + ".unclassified_sequences_2.fastq"
                classified_output_1 = kraken_out_prefix + ".classified_sequences_1.fastq"
                classified_output_2 = kraken_out_prefix + ".classified_sequences_2.fastq"
                kraken_outputs.extend((unclassified_output_1, unclassified_output_2, classified_output_1, classified_output_2))

            elif readset.run_type == "SINGLE_END":
                unclassified_output = kraken_out_prefix + ".unclassified_sequences.fastq"
                classified_output = kraken_out_prefix + ".classified_sequences.fastq"
                kraken_outputs.extend((unclassified_output, classified_output))

        quast_outputs = []
        flagstat_outputs = []
        bedgraph_outputs = []
        picard_outputs = []
        for sample in self.samples:
            # Finding quast outputs
            ivar_quast_output_dir = os.path.join("metrics", "dna", sample.name, "quast_metrics_ivar")
            quast_outputs.extend([
                os.path.join(ivar_quast_output_dir, "report.html"),
                os.path.join(ivar_quast_output_dir, "report.pdf"),
                os.path.join(ivar_quast_output_dir, "report.tex"),
                os.path.join(ivar_quast_output_dir, "report.tsv"),
                os.path.join(ivar_quast_output_dir, "report.txt")
                ])
            # Finding flagstat outputs
            flagstat_directory = os.path.join("metrics", "dna", sample.name, "flagstat")
            alignment_directory = os.path.join("alignment", sample.name)
            input_bams = [
                os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.bam"),
                os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam"),
                os.path.join(alignment_directory, sample.name + ".sorted.bam")
            ]
            for input_bam in input_bams:
                flagstat_outputs.append(os.path.join(flagstat_directory, re.sub("\.bam$", ".flagstat", os.path.basename(input_bam))))
            # Finding BedGraph file
            [input_bam] = self.select_input_files([
                [os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.bam")]
            ])
            bedgraph_outputs.append(os.path.join(alignment_directory, re.sub("\.bam$", ".BedGraph", os.path.basename(input_bam))))

            # Picard collect multiple metrics outputs
            picard_directory = os.path.join("metrics", "dna", sample.name, "picard_metrics")
            picard_out = os.path.join(picard_directory, re.sub("\.bam$", "", os.path.basename(input_bam)) + ".all.metrics")

            if library[sample] == "PAIRED_END":
                picard_outputs.extend([
                    picard_out + ".alignment_summary_metrics",
                    picard_out + ".insert_size_metrics"
                ])
            else:
                picard_outputs.extend([
                    picard_out + ".alignment_summary_metrics",
                ])

        covid_collect_metrics_inputs = []
        covid_collect_metrics_inputs.extend(kraken_outputs)
        covid_collect_metrics_inputs.append(input_bam)
        covid_collect_metrics_inputs.extend(quast_outputs)
        covid_collect_metrics_inputs.extend(flagstat_outputs)
        covid_collect_metrics_inputs.extend(bedgraph_outputs)
        covid_collect_metrics_inputs.extend(picard_outputs)
        jobs.append(
            concat_jobs([
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
                    command="""\\
echo "Preparing to run metadata..." && \\
echo "run_name,{run_name}
genpipes_version,{genpipes_version}
cluster_server,{cluster_server}
assembly_synonyms,{assembly_synonyms}
sequencing_technology,{sequencing_technology}" > {run_metadata} && \\
echo "Software Versions
{modules_all}" > {software_version}""".format(
    output_dir=self.output_dir,
    run_name=config.param('prepare_report', 'run_name', required=True),
    genpipes_version=self.genpipes_version.strip(),
    cluster_server=config.param('prepare_report', 'cluster_server'),
    assembly_synonyms=config.param('prepare_report', 'assembly_synonyms'),
    sequencing_technology=config.param('prepare_report', 'sequencing_technology'),
    run_metadata=run_metadata,
    modules_all="\n".join(modules),
    software_version=software_version
    )
                    )
                ],
                name="prepare_table." + config.param('prepare_report', 'run_name', required=True)
                )
            )

        return jobs

    def prepare_report_ivar(self):

        jobs = []

        readset_file=os.path.relpath(self.args.readsets.name, self.output_dir)
        ivar_readset_file_report="report.readset_ivar.tsv"

        software_version = os.path.join("report", "software_versions.csv")
        run_metadata = os.path.join("report", "run_metadata.csv")

        ivar_ncovtools_directory = os.path.join("report", "ncov_tools_ivar")
        ivar_metadata = os.path.join(ivar_ncovtools_directory, "metadata.tsv")
        ivar_ncovtools_data_directory = os.path.join(ivar_ncovtools_directory, "data")
        ivar_ncovtools_config = os.path.join(ivar_ncovtools_directory, "config.yaml")

        modules = []
        # Retrieve all unique module version values in config files
        # assuming that all module key names start with "module_"
        for section in config.sections():
            for name, value in config.items(section):
                if re.search("^module_", name) and value not in modules:
                    modules.append(value)

        job = concat_jobs([
            bash.mkdir(ivar_ncovtools_data_directory),
            Job(
                    input_files=[],
                    output_files=[ivar_readset_file_report, ivar_metadata],
                    command="""\\
head -n 1 {readset_file} > {ivar_readset_file_report} && \\
echo -e "sample\\tct\\tdate" > {ivar_metadata}""".format(
    readset_file=readset_file,
    ivar_readset_file_report=ivar_readset_file_report,
    ivar_metadata=ivar_metadata
    )
                ),
            ])

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            filtered_bam = os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")
            primer_trimmed_bam = os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.bam")
            consensus_directory = os.path.join("consensus", sample.name)
            ivar_consensus = os.path.join(consensus_directory, sample.name + ".consensus.fasta")
            variant_directory = os.path.join("variant", sample.name)
            ivar_variants = os.path.join(variant_directory, sample.name + ".variants.tsv")

            ivar_output_filtered_bam = os.path.join(ivar_ncovtools_data_directory, os.path.basename(filtered_bam))
            ivar_output_primer_trimmed_bam = os.path.join(ivar_ncovtools_data_directory, os.path.basename(primer_trimmed_bam))
            ivar_output_consensus = os.path.join(ivar_ncovtools_data_directory, os.path.basename(ivar_consensus))
            ivar_output_variants = os.path.join(ivar_ncovtools_data_directory, os.path.basename(ivar_variants))

            job = concat_jobs([
                        job,
                        Job(
                            input_files=[filtered_bam, primer_trimmed_bam, ivar_consensus, ivar_variants],
                            output_files=[ivar_output_filtered_bam, ivar_output_primer_trimmed_bam, ivar_output_consensus, ivar_output_variants],
                            command="""\\
echo "Linking files for ncov_tools for sample {sample_name}..." && \\
if [ "$(ls -1 {filtered_bam})" != "" ] && [ "$(ls -1 {primer_trimmed_bam})" != "" ] && [ "$(ls -1 {ivar_consensus})" != "" ] && [ "$(ls -1 {ivar_variants})" != "" ];
  then
    ln -fs $(pwd -P )/$(ls -1 {filtered_bam}) {ivar_output_filtered_bam} && \\
    ln -fs $(pwd -P )/$(ls -1 {primer_trimmed_bam}) {ivar_output_primer_trimmed_bam} && \\
    ln -fs $(pwd -P )/$(ls -1 {ivar_consensus}) {ivar_output_consensus} && \\
    ln -fs $(pwd -P )/$(ls -1 {ivar_variants}) {ivar_output_variants} && \\
    grep {sample_name} {readset_file} >> {ivar_readset_file_report} && \\
    echo -e "{sample_name}\\tNA\\tNA" >> {ivar_metadata}
fi""".format(
    readset_file=readset_file,
    ivar_readset_file_report=ivar_readset_file_report,
    filtered_bam=filtered_bam,
    primer_trimmed_bam=primer_trimmed_bam,
    ivar_consensus=ivar_consensus,
    ivar_variants=ivar_variants,
    ivar_output_filtered_bam=ivar_output_filtered_bam,
    ivar_output_primer_trimmed_bam=ivar_output_primer_trimmed_bam,
    ivar_output_consensus=ivar_output_consensus,
    ivar_output_variants=ivar_output_variants,
    sample_name=sample.name,
    ivar_metadata=ivar_metadata
    )
                            )
                    ],
                    samples=[sample]
                    )
        jobs.append(
            concat_jobs([
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
                Job(input_files=[],
                    output_files=[],
                    module_entries=[
                        ['prepare_report', 'module_R'],
                        ['prepare_report', 'module_CoVSeQ_tools']
                    ],
                    command="""\\
module purge && \\
module load {R_covseqtools}""".format(
    R_covseqtools=config.param('prepare_report', 'module_R') + " " + config.param('prepare_report', 'module_CoVSeQ_tools'),
    output_dir=self.output_dir)
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
                name="prepare_report." + config.param('prepare_report', 'run_name', required=True)
                )
            )

        return jobs

    def prepare_report_freebayes(self):

        jobs = []

        readset_file=os.path.relpath(self.args.readsets.name, self.output_dir)
        freebayes_readset_file_report="report.readset_freebayes.tsv"

        software_version = os.path.join("report", "software_versions.csv")
        run_metadata = os.path.join("report", "run_metadata.csv")

        freebayes_ncovtools_directory = os.path.join("report", "ncov_tools_freebayes")
        freebayes_metadata = os.path.join(freebayes_ncovtools_directory, "metadata.tsv")
        freebayes_ncovtools_data_directory = os.path.join(freebayes_ncovtools_directory, "data")
        freebayes_ncovtools_config = os.path.join(freebayes_ncovtools_directory, "config.yaml")

        job = concat_jobs([
            bash.mkdir(freebayes_ncovtools_data_directory),
            Job(
                    input_files=[],
                    output_files=[freebayes_readset_file_report, freebayes_metadata],
                    command="""\\
head -n 1 {readset_file} > {freebayes_readset_file_report} && \\
echo -e "sample\\tct\\tdate" > {freebayes_metadata}""".format(
    readset_file=readset_file,
    freebayes_readset_file_report=freebayes_readset_file_report,
    freebayes_metadata=freebayes_metadata
    )
                ),
            ])

        for sample in self.samples:
        
            filtered_bam = os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")
            primer_trimmed_bam = os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.bam")
            consensus_directory = os.path.join("consensus", sample.name)
            freebayes_consensus = os.path.join(consensus_directory, sample.name + ".freebayes_calling.consensus.renamed.fasta")
            variant_directory = os.path.join("variant", sample.name)
            freebayes_variants = os.path.join(variant_directory, sample.name + ".freebayes_calling.consensus.vcf")

            freebayes_output_filtered_bam = os.path.join(freebayes_ncovtools_data_directory, os.path.basename(filtered_bam))
            freebayes_output_primer_trimmed_bam = os.path.join(freebayes_ncovtools_data_directory, os.path.basename(primer_trimmed_bam))
            freebayes_output_consensus = os.path.join(freebayes_ncovtools_data_directory, os.path.basename(freebayes_consensus))
            freebayes_output_variants = os.path.join(freebayes_ncovtools_data_directory, os.path.basename(freebayes_variants))

            job = concat_jobs([
                        job,
                        Job(
                            input_files=[filtered_bam, primer_trimmed_bam, freebayes_consensus, freebayes_variants],
                            output_files=[freebayes_output_filtered_bam, freebayes_output_primer_trimmed_bam, freebayes_output_consensus, freebayes_output_variants],
                            command="""\\
echo "Linking files for ncov_tools for sample {sample_name}..." && \\
if [ "$(ls -1 {filtered_bam})" != "" ] && [ "$(ls -1 {primer_trimmed_bam})" != "" ] && [ "$(ls -1 {freebayes_consensus})" != "" ] && [ "$(ls -1 {freebayes_variants})" != "" ];
  then
    ln -fs $(pwd -P )/$(ls -1 {filtered_bam}) {freebayes_output_filtered_bam} && \\
    ln -fs $(pwd -P )/$(ls -1 {primer_trimmed_bam}) {freebayes_output_primer_trimmed_bam} && \\
    ln -fs $(pwd -P )/$(ls -1 {freebayes_consensus}) {freebayes_output_consensus} && \\
    ln -fs $(pwd -P )/$(ls -1 {freebayes_variants}) {freebayes_output_variants} && \\
    grep {sample_name} {readset_file} >> {freebayes_readset_file_report} && \\
    echo -e "{sample_name}\\tNA\\tNA" >> {freebayes_metadata}
fi""".format(
    readset_file=readset_file,
    freebayes_readset_file_report=freebayes_readset_file_report,
    filtered_bam=filtered_bam,
    primer_trimmed_bam=primer_trimmed_bam,
    sample_name=sample.name,
    freebayes_consensus=freebayes_consensus,
    freebayes_variants=freebayes_variants,
    freebayes_output_filtered_bam=freebayes_output_filtered_bam,
    freebayes_output_primer_trimmed_bam=freebayes_output_primer_trimmed_bam,
    freebayes_output_consensus=freebayes_output_consensus,
    freebayes_output_variants=freebayes_output_variants,
    freebayes_metadata=freebayes_metadata
    )
                            )
                    ],
                    samples=[sample]
                    )

        jobs.append(
            concat_jobs([
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
                Job(input_files=[],
                    output_files=[],
                    module_entries=[
                        ['prepare_report', 'module_R'],
                        ['prepare_report', 'module_CoVSeQ_tools']
                    ],
                    command="""\\
module purge && \\
module load {R_covseqtools}""".format(
    R_covseqtools=config.param('prepare_report', 'module_R') + " " + config.param('prepare_report', 'module_CoVSeQ_tools'),
    output_dir=self.output_dir)
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
                name="prepare_report." + config.param('prepare_report', 'run_name', required=True)
                )
            )

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
            self.rename_consensus_header_ivar,
            self.rename_consensus_header_freebayes,
            self.ncovtools_quickalign,
            self.prepare_table,
            self.prepare_report_ivar,
            self.prepare_report_freebayes
            # self.run_multiqc
        ]

if __name__ == '__main__':
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        CoVSeQ()
