#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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
import argparse
import csv
import logging
import math
import os
import re
import subprocess
import string
import sys
import time

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs
import utils.utils

from pipelines import common

from bfx.sequence_dictionary import parse_sequence_dictionary_file, split_by_size

from bfx import bvatools
from bfx import bwa
from bfx import gatk4
from bfx import gq_seq_utils
from bfx import homer
from bfx import macs2
from bfx import multiqc
from bfx import picard
from bfx import sambamba
from bfx import samtools
from bfx import tools
from bfx import trimmomatic
from bfx import ucsc
from bfx import differential_binding
# from pipelines.dnaseq import dnaseq

from bfx import bash_cmd as bash

from bfx.readset import parse_illumina_readset_file
from bfx.design import parse_chipseq_design_file

log = logging.getLogger(__name__)


class ChipSeq(common.Illumina):
    """
    ChIP-Seq Pipeline
    =================

    ChIP-Seq experiments allows the Isolation and sequencing of genomic DNA bound by a specific transcription factor,
    covalently modified histone, or other nuclear protein. The pipeline starts by trimming adaptors and
    low quality bases and mapping the reads (single end or paired end ) to a reference genome using bwa.
    Reads are filtered by mapping quality and duplicate reads are marked. Then, Homer quality control routines
    are used to provide information and feedback about the quality of the experiment. Peak calls is executed by MACS
    and annotation and motif discovery for narrow peaks are executed using Homer. Statistics of annotated peaks
    are produced for narrow peaks and a standard report is generated.

    An example of the ChIP-Seq report for an analysis on public ENCODE data is available for illustration purpose only:
    [ChIP-Seq report](http://gqinnovationcenter.com/services/bioinformatics/tools/chipReport/index.html).

    [Here](https://bitbucket.org/mugqic/mugqic_pipelines/downloads/MUGQIC_Bioinfo_ChIP-Seq.pptx)
    is more information about ChIP-Seq pipeline that you may find interesting.
    """

    def __init__(self, protocol="chipseq"):
        self._protocol = protocol
        # Add pipeline specific arguments
        self.argparser.add_argument("-d", "--design", help="design file", type=argparse.FileType('r'), required=False)
        self.argparser.add_argument("-t", "--type", help="Type of pipeline (default chipseq)", choices=["chipseq", "atacseq"], default="chipseq")
        super(ChipSeq, self).__init__(protocol)

    @property
    def output_dirs(self):
        dirs = {'alignment_output_directory': 'alignment',
                'report_output_directory': 'report',
                'metrics_output_directory': 'metrics',
                'homer_output_directory': 'tags',
                'graphs_output_directory': 'graphs',
                'tracks_output_directory': 'tracks',
                'macs_output_directory': 'peak_call',
                'anno_output_directory': 'annotation',
                'ihecA_output_directory': 'ihec_alignment',
                'ihecM_output_directory': 'ihec_metrics',
                'dba_output_directory': 'differential_binding'
                }
        return dirs

    @property
    def mark_type_conversion(self):
        dirs = {'N': 'narrow',
                'B': 'broad',
                'I': 'Input'
                }
        return dirs

    @property
    def ucsc_genome(self):
        genome_source = config.param('DEFAULT', 'source')
        if genome_source == "UCSC":
            genome = config.param('DEFAULT', 'assembly')
        else:
            genome = config.param('DEFAULT', 'assembly_synonyms')
        return genome

    @property
    def readsets(self):
        flag = False
        if not hasattr(self, "_readsets"):
            if self.args.readsets:
                self._readsets = parse_illumina_readset_file(self.args.readsets.name)
                for readset in self.readsets:
                    if not readset.mark_name:
                        _raise(SanitycheckError("Error: missing readset MarkName for " + readset.name))
                        flag = True
                    elif not readset.mark_type:
                        _raise(SanitycheckError("Error: missing readset MarkType for " + readset.name))
                        flag = True
                if flag:
                    exit()
            else:
                self.argparser.error("argument -r/--readsets is required!")

        return self._readsets

    @property
    def contrasts(self):
        flag = False

        if self.args.design:
            self._contrast = parse_chipseq_design_file(self.args.design.name, self.samples)
        else:
            self.argparser.error("argument -d/--design is required!")

        return self._contrast

    def sequence_dictionary_variant(self):
        if not hasattr(self, "_sequence_dictionary_variant"):
            self._sequence_dictionary_variant = parse_sequence_dictionary_file(config.param('DEFAULT', 'genome_dictionary', param_type='filepath'), variant=True)
        return self._sequence_dictionary_variant

    def mappable_genome_size(self):
        genome_index = csv.reader(open(config.param('DEFAULT', 'genome_fasta', param_type='filepath') + ".fai", 'r'), delimiter='\t')
        # 2nd column of genome index contains chromosome length
        # HOMER and MACS2 mappable genome size (without repetitive features) is about 80 % of total size
        return int(sum([int(chromosome[1]) for chromosome in genome_index]) * config.param('DEFAULT', 'mappable_genome_size', param_type='float', required=True))

    def trimmomatic(self):
        """
        Raw reads quality trimming and removing of Illumina adapters is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
        If an adapter FASTA file is specified in the config file (section 'trimmomatic', param 'adapter_fasta'),
        it is used first. Else, 'Adapter1' and 'Adapter2' columns from the readset file are used to create
        an adapter FASTA file, given then to Trimmomatic. For PAIRED_END readsets, readset adapters are
        reversed-complemented and swapped, to match Trimmomatic Palindrome strategy. For SINGLE_END readsets,
        only Adapter1 is used and left unchanged.

        This step takes as input files:

        1. FASTQ files from the readset file if available
        2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """
        jobs = []
        for readset in self.readsets:
            # log.info(readset.mark_name)
            trim_directory = os.path.join("trim", readset.sample.name, readset.mark_name)
            trim_file_prefix = os.path.join(trim_directory, readset.name + ".trim.")
            trim_log = trim_file_prefix + "log"

            # Use adapter FASTA in config file if any, else create it from readset file
            adapter_fasta = config.param('trimmomatic', 'adapter_fasta', required=False, param_type='filepath')
            adapter_job = None
            if not adapter_fasta:
                adapter_fasta = trim_file_prefix + "adapters.fa"
                if readset.run_type == "PAIRED_END":
                    if readset.adapter1 and readset.adapter2:
                        # WARNING: Reverse-complement and swap readset adapters for Trimmomatic Palindrome strategy
                        adapter_job = Job(command="""\
`cat > {adapter_fasta} << END
>Prefix/1
{sequence1}
>Prefix/2
{sequence2}
END
`""".format(adapter_fasta=adapter_fasta,
            sequence1=readset.adapter2.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1],
            sequence2=readset.adapter1.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]))
                    else:
                        _raise(SanitycheckError(
                            "Error: missing adapter1 and/or adapter2 for PAIRED_END readset \"" + readset.name + "\", or missing adapter_fasta parameter in config file!"))
                elif readset.run_type == "SINGLE_END":
                    if readset.adapter1:
                        adapter_job = Job(command="""\
`cat > {adapter_fasta} << END
>Single
{sequence}
END
`""".format(adapter_fasta=adapter_fasta, sequence=readset.adapter1))
                    else:
                        _raise(SanitycheckError(
                            "Error: missing adapter1 for SINGLE_END readset \"" + readset.name + "\", or missing adapter_fasta parameter in config file!"))

            trim_stats = trim_file_prefix + "stats.csv"
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[readset.fastq1, readset.fastq2]]
                if readset.bam:
                    candidate_fastq1 = os.path.join(self.output_dir, "raw_reads", readset.sample.name,
                                                    readset.name + ".pair1.fastq.gz")
                    candidate_fastq2 = os.path.join(self.output_dir, "raw_reads", readset.sample.name,
                                                    readset.name + ".pair2.fastq.gz")
                    candidate_input_files.append([candidate_fastq1, candidate_fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                job = trimmomatic.trimmomatic(
                    fastq1,
                    fastq2,
                    trim_file_prefix + "pair1.fastq.gz",
                    trim_file_prefix + "single1.fastq.gz",
                    trim_file_prefix + "pair2.fastq.gz",
                    trim_file_prefix + "single2.fastq.gz",
                    None,
                    readset.quality_offset,
                    adapter_fasta,
                    trim_log
                )
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[readset.fastq1]]
                if readset.bam:
                    candidate_input_files.append([os.path.join(self.output_dir, "raw_reads", readset.sample.name,
                                                               readset.name + ".single.fastq.gz")])
                [fastq1] = self.select_input_files(candidate_input_files)
                job = trimmomatic.trimmomatic(
                    fastq1,
                    None,
                    None,
                    None,
                    None,
                    None,
                    trim_file_prefix + "single.fastq.gz",
                    readset.quality_offset,
                    adapter_fasta,
                    trim_log
                )
            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                                        "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            if adapter_job:
                job = concat_jobs([adapter_job, job])

            jobs.append(concat_jobs([
                # Trimmomatic does not create output directory by default
                bash.mkdir(trim_directory),
                job
            ], name="trimmomatic." + readset.name, samples=[readset.sample]))
        return jobs

    def merge_trimmomatic_stats(self):
        """
        The trim statistics per readset are merged at this step.
        """

        read_type = "Paired" if self.run_type == 'PAIRED_END' else "Single"
        readset_merge_trim_stats = os.path.join("metrics", "trimReadsetTable.tsv")
        job = concat_jobs([
            bash.mkdir(self.output_dirs['metrics_output_directory']),
            Job(
                command="""
echo -e "Sample\\tReadset\\tMark Name\\tRaw {read_type} Reads #\\tSurviving {read_type} Reads #\\tSurviving {read_type} Reads %" > {readset_merge_trim_stats}""".format(
                    read_type=read_type,
                    readset_merge_trim_stats=readset_merge_trim_stats
                )
            )
        ])

        for readset in self.readsets:
            trim_log = os.path.join("trim", readset.sample.name, readset.mark_name, readset.name + ".trim.log")
            if readset.run_type == "PAIRED_END":
                # Retrieve readset raw and surviving reads from trimmomatic log using ugly Perl regexp
                perl_command = "perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/{readset.sample.name}\\t{readset.name}\\t{readset.mark_name}\\t\\1\\t\\2/'".format(
                    readset=readset)
            elif readset.run_type == "SINGLE_END":
                perl_command = "perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/{readset.sample.name}\\t{readset.name}\\t{readset.mark_name}\\t\\1\\t\\2/'".format(
                    readset=readset)

            job = concat_jobs([
                job,
                Job(
                    [trim_log],
                    [readset_merge_trim_stats],
                    module_entries=[['merge_trimmomatic_stats', 'module_perl']],
                    # Create readset trimming stats TSV file with paired or single read count using ugly awk
                    command="""\
grep ^Input {trim_log} | \\
{perl_command} | \\
awk '{{OFS="\\t"; print $0, $5 / $4 * 100}}' \\
  >> {readset_merge_trim_stats}""".format(
                        trim_log=trim_log,
                        perl_command=perl_command,
                        readset_merge_trim_stats=readset_merge_trim_stats
                    ),
                    samples=[readset.sample]
                )
            ])

        sample_merge_trim_stats = os.path.join("metrics", "trimSampleTable.tsv")
        report_file = os.path.join("report", "Illumina.merge_trimmomatic_stats.md")
        return [concat_jobs([
            job,
            Job(
                [readset_merge_trim_stats],
                [sample_merge_trim_stats],
                # Create sample trimming stats TSV file with total read counts (i.e. paired * 2 if applicable) using ugly awk
                command="""\
cut -f1,3- {readset_merge_trim_stats} | awk -F"\\t" '{{OFS="\\t"; if (NR==1) {{if ($3=="Raw Paired Reads #") {{paired=1}};print "Sample", "Mark Name", "Raw Reads #", "Surviving Reads #", "Surviving %"}} else {{if (paired) {{$3=$3*2; $4=$4*2}}; sample[$1$2]=$1; markname[$1$2]=$2; raw[$1$2]+=$3; surviving[$1$2]+=$4}}}}END{{for (samplemark in raw){{print sample[samplemark], markname[samplemark], raw[samplemark], surviving[samplemark], surviving[samplemark] / raw[samplemark] * 100}}}}' \\
  > {sample_merge_trim_stats}""".format(
                    readset_merge_trim_stats=readset_merge_trim_stats,
                    sample_merge_trim_stats=sample_merge_trim_stats
                )
            ),
            Job(
                [sample_merge_trim_stats],
                [report_file],
                [['merge_trimmomatic_stats', 'module_pandoc']],
                command="""\
mkdir -p report && \\
cp {readset_merge_trim_stats} {sample_merge_trim_stats} report/ && \\
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "\\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----|-----|-----:|-----:|-----:"}} else {{print $1, $2, $3, sprintf("%\\47d", $4), sprintf("%\\47d", $5), sprintf("%.1f", $6)}}}}' {readset_merge_trim_stats}` && \\
pandoc \\
  {report_template_dir}/{basename_report_file} \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable trailing_min_quality={trailing_min_quality} \\
  --variable min_length={min_length} \\
  --variable read_type={read_type} \\
  --variable trim_readset_table="$trim_readset_table_md" \\
  --to markdown \\
  > {report_file}""".format(
                    trailing_min_quality=config.param('trimmomatic', 'trailing_min_quality', param_type='int'),
                    min_length=config.param('trimmomatic', 'min_length', param_type='posint'),
                    read_type=read_type,
                    report_template_dir=self.report_template_dir,
                    readset_merge_trim_stats=readset_merge_trim_stats,
                    sample_merge_trim_stats=sample_merge_trim_stats,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file]
                )
            ],
            name="merge_trimmomatic_stats."
            )
        ]
        # TODO: replace ".".join([sample.name for sample in self.samples])) by timestamp to avoid too long naming issue

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
            trim_directory = os.path.join("trim", readset.sample.name, readset.mark_name)
            trim_file_prefix = os.path.join(trim_directory, readset.name)
            alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], readset.sample.name,
                                               readset.mark_name)
            readset_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam")
            index_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam.bai")

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [
                    [trim_file_prefix + ".trim.pair1.fastq.gz", trim_file_prefix + ".trim.pair2.fastq.gz"]
                ]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [
                    [trim_file_prefix + ".trim.single.fastq.gz"]
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
                                       (
                                           "\\tPU:run" + readset.run + "_" + readset.lane if readset.run and readset.lane else "") + \
                                       ("\\tCN:" + config.param('mapping_bwa_mem_sambamba',
                                                                'sequencing_center') if config.param(
                                           'mapping_bwa_mem_sambamba', 'sequencing_center', required=False) else "") + \
                                       ("\\tPL:" + config.param('mapping_bwa_mem_sambamba',
                                                                'sequencing_technology') if config.param(
                                           'mapping_bwa_mem_sambamba', 'sequencing_technology',
                                           required=False) else "Illumina") + \
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
                            other_options=config.param('mapping_bwa_mem_sambamba', 'sambamba_sort_other_options',
                                                       required=False)
                        )
                    ]),
                    sambamba.index(
                        readset_bam,
                        index_bam,
                        other_options=config.param('mapping_bwa_mem_sambamba', 'sambamba_index_other_options',
                                                   required=False)
                    )
                ],
                    name="mapping_bwa_mem_sambamba." + readset.name,
                    samples=[readset.sample]
                )
            )

        return jobs


    def sambamba_merge_bam_files(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [Sambamba]().

        This step takes as input files:

        1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
        2. Else, BAM files from the readset file
        """

        jobs = []
        for sample in self.samples:
            for mark_name in sample.marks:
                alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], sample.name,
                                                   mark_name)
                # Find input readset BAMs first from previous bwa_mem_picard_sort_sam job, then from original BAMs in the readset sheet.
                readset_bams = [os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam") for
                                readset in sample.readsets if readset.mark_name == mark_name]
                sample_bam = os.path.join(alignment_directory, sample.name + "." + mark_name + ".sorted.bam")

                # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
                if len(readset_bams) == 1:
                    readset_bam = readset_bams[0]
                    readset_index = re.sub("\.bam$", ".bam.bai", readset_bam)
                    sample_index = re.sub("\.bam$", ".bam.bai", sample_bam)

                    jobs.append(
                        concat_jobs([
                            bash.mkdir(
                                os.path.dirname(sample_bam)
                            ),
                            bash.ln(
                                readset_bam,
                                sample_bam,
                                self.output_dir
                            ),
                            bash.ln(
                                readset_index,
                                sample_index,
                                self.output_dir
                            )
                        ],
                            name="symlink_readset_sample_bam." + sample.name + "." + mark_name,
                            samples=[sample]
                        )
                    )

                elif len(sample.readsets) > 1:
                    jobs.append(
                        concat_jobs([
                            bash.mkdir(
                                os.path.dirname(sample_bam)
                            ),
                            sambamba.merge(
                                readset_bams,
                                sample_bam,
                                ini_section="sambamba_merge_bam_files"
                            )
                        ],
                            name="sambamba_merge_bam_files." + sample.name + "." + mark_name,
                            samples=[sample]
                        )
                    )
        return jobs


    def sambamba_mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using [Sambamba]().
        """

        jobs = []
        for sample in self.samples:
            for mark_name in sample.marks:
                alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], sample.name,
                                                   mark_name)
                input_bam = os.path.join(alignment_directory, sample.name + "." + mark_name + ".sorted.bam")
                output_bam = os.path.join(alignment_directory, sample.name + "." + mark_name + ".sorted.dup.bam")
                # metrics_file = alignment_file_prefix + ".sorted.dup.metrics"

                jobs.append(
                    concat_jobs([
                        bash.mkdir(os.path.dirname(output_bam)),
                        sambamba.markdup(
                            input_bam,
                            output_bam,
                            tmp_dir=config.param('sambamba_mark_duplicates', 'tmp_dir', required=True),
                            other_options=config.param('sambamba_mark_duplicates', 'other_options', required=False)
                        )
                    ],
                        name="sambamba_mark_duplicates." + sample.name + "." + mark_name,
                        samples=[sample]
                    )
                )

        report_file = os.path.join(self.output_dirs['report_output_directory'], "ChipSeq.sambamba_mark_duplicates.md")
        jobs.append(
            Job(
                [os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name,
                              sample.name + "." + mark_name + ".sorted.dup.bam") for sample in self.samples for
                 mark_name in sample.marks],
                [report_file],
                command="""\
mkdir -p {report_dir} && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file,
                    report_dir=self.output_dirs['report_output_directory']
                ),
                report_files=[report_file],
                name="sambamba_mark_duplicates_report"#".".join([sample.name for sample in self.samples])
                )
        )

        return jobs

    def sambamba_view_filter(self):
        """
        Filter out unmapped reads and low quality reads [Sambamba](http://www.htslib.org/).
        """

        jobs = []
        for sample in self.samples:
            for mark_name in sample.marks:
                alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], sample.name,
                                                   mark_name)
                input_bam = os.path.join(alignment_directory, sample.name + "." + mark_name + ".sorted.dup.bam")
                output_bam = os.path.join(alignment_directory,
                                          sample.name + "." + mark_name + ".sorted.dup.filtered.bam")
                output_bam_index = os.path.join(alignment_directory,
                                                sample.name + "." + mark_name + ".sorted.dup.filtered.bam.bai")
                jobs.append(
                    concat_jobs([
                        bash.mkdir(os.path.dirname(output_bam)),
                        sambamba.view(
                            input_bam,
                            output_bam,
                            """-t {threads} -f bam -F \"not unmapped and not failed_quality_control and mapping_quality >= {min_mapq}\"""".format(
                                threads=config.param('sambamba_view_filter', 'threads'),
                                min_mapq=config.param('sambamba_view_filter', 'min_mapq'))
                        ),
                        sambamba.index(
                            output_bam,
                            output_bam_index
                        )
                    ],
                        name="sambamba_view_filter." + sample.name + "." + mark_name,
                        samples=[sample]
                    )
                )
        report_file = os.path.join(self.output_dirs['report_output_directory'], "ChipSeq.sambamba_view_filter.md")
        jobs.append(
            Job(
                [os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name,
                              sample.name + "." + mark_name + ".sorted.dup.filtered.bam") for sample in self.samples for
                 mark_name in sample.marks],
                [report_file],
                [['sambamba_view_filter', 'module_pandoc']],
                command="""\
mkdir -p {report_dir} && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable min_mapq="{min_mapq}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    min_mapq=config.param('sambamba_view_filter', 'min_mapq', param_type='int'),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file,
                    report_dir=self.output_dirs['report_output_directory']
                ),
                report_files=[report_file],
                name="sambamba_view_filter_report"#".".join([sample.name for sample in self.samples])
                )
        )

        return jobs


    def metrics(self):
        """
        The number of raw/filtered and aligned reads per sample are computed at this stage.
        """

        jobs = []

        samples_associative_array = []

        metrics_output_directory = self.output_dirs['metrics_output_directory']

        inputs_report = []

        for sample in self.samples:
            samples_associative_array.append("[\"" + sample.name + "\"]=\"" + " ".join(sample.marks.keys()) + "\"")
            for mark_name in sample.marks:
                alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], sample.name,
                                                   mark_name)
                raw_bam_file = os.path.join(alignment_directory, sample.name + "." + mark_name + ".sorted.dup.bam")
                bam_file = os.path.join(alignment_directory,
                                        sample.name + "." + mark_name + ".sorted.dup.filtered.bam")

                jobs.append(
                    concat_jobs([
                        bash.mkdir(os.path.join(metrics_output_directory, sample.name, mark_name)),
                        picard.collect_multiple_metrics(
                            bam_file,
                            os.path.join(metrics_output_directory, sample.name, mark_name,
                                         re.sub("bam$", "all.metrics", os.path.basename(bam_file))),
                            library_type=self.run_type
                        )
                    ],
                        name="picard_collect_multiple_metrics." + sample.name + "." + mark_name
                    )
                )

                jobs.append(
                    concat_jobs([
                        bash.mkdir(os.path.join(metrics_output_directory, sample.name, mark_name)),
                        sambamba.flagstat(
                            raw_bam_file,
                            os.path.join(metrics_output_directory, sample.name, mark_name,
                                         re.sub("\.bam$", ".flagstat", os.path.basename(raw_bam_file)))
                        ),
                        sambamba.flagstat(
                            bam_file,
                            os.path.join(metrics_output_directory, sample.name, mark_name,
                                         re.sub("\.bam$", ".flagstat", os.path.basename(bam_file)))
                        )
                    ],
                        name="metrics_flagstat." + sample.name + "." + mark_name
                    )
                )
                inputs_report.extend((os.path.join(metrics_output_directory, sample.name, mark_name,
                                                   re.sub("\.bam$", ".flagstat", os.path.basename(raw_bam_file))),
                                      os.path.join(metrics_output_directory, sample.name, mark_name,
                                                   re.sub("\.bam$", ".flagstat", os.path.basename(bam_file))),
                                      os.path.join(self.output_dirs['alignment_output_directory'], sample.name,
                                                   mark_name,
                                                   sample.name + "." + mark_name + ".sorted.dup.filtered.bam")))

        trim_metrics_file = os.path.join(metrics_output_directory, "trimSampleTable.tsv")
        metrics_file = os.path.join(metrics_output_directory, "SampleMetrics.tsv")
        report_metrics_file = os.path.join(self.output_dirs['report_output_directory'], "SampleMetrics.tsv")
        report_file = os.path.join(self.output_dirs['report_output_directory'], "ChipSeq.metrics.md")
        jobs.append(
            Job(
                inputs_report,
                [report_metrics_file],
                [['metrics', 'module_pandoc']],
                # Retrieve number of aligned and duplicate reads from sample flagstat files
                # Merge trimming stats per sample with aligned and duplicate stats using ugly awk
                # Format merge stats into markdown table using ugly awk (knitr may do this better)
                command="""\
module load {sambamba} && \\
mkdir -p {metrics_dir}
cp /dev/null {metrics_file} && \\
declare -A samples_associative_array=({samples_associative_array}) && \\
for sample in ${{!samples_associative_array[@]}}
do
  for mark_name in ${{samples_associative_array[$sample]}}
  do
    raw_flagstat_file={metrics_dir}/$sample/$mark_name/$sample.$mark_name.sorted.dup.flagstat
    filtered_flagstat_file={metrics_dir}/$sample/$mark_name/$sample.$mark_name.sorted.dup.filtered.flagstat
    bam_file={alignment_dir}/$sample/$mark_name/$sample.$mark_name.sorted.dup.filtered.bam
    raw_supplementarysecondary_reads=`bc <<< $(grep "secondary" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* secondary.*//')+$(grep "supplementary" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* supplementary.*//')`
    mapped_reads=`bc <<< $(grep "mapped (" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* mapped (.*)//')-$raw_supplementarysecondary_reads`
    filtered_supplementarysecondary_reads=`bc <<< $(grep "secondary" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* secondary.*//')+$(grep "supplementary" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* supplementary.*//')`
    filtered_reads=`bc <<< $(grep "in total" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$filtered_supplementarysecondary_reads`
    filtered_mapped_reads=`bc <<< $(grep "mapped (" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* mapped (.*)//')-$filtered_supplementarysecondary_reads`
    filtered_mapped_rate=`echo "scale=4; 100*$filtered_mapped_reads/$filtered_reads" | bc -l`
    filtered_dup_reads=`grep "duplicates" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* duplicates$//'`
    filtered_dup_rate=`echo "scale=4; 100*$filtered_dup_reads/$filtered_mapped_reads" | bc -l`
    filtered_dedup_reads=`echo "$filtered_mapped_reads-$filtered_dup_reads" | bc -l`
    if [[ -s {trim_metrics_file} ]]
      then
        raw_reads=$(grep -P "${{sample}}\\t${{mark_name}}" {trim_metrics_file} | cut -f 3)
        raw_trimmed_reads=`bc <<< $(grep "in total" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$raw_supplementarysecondary_reads`
        mapped_reads_rate=`echo "scale=4; 100*$mapped_reads/$raw_trimmed_reads" | bc -l`
        raw_trimmed_rate=`echo "scale=4; 100*$raw_trimmed_reads/$raw_reads" | bc -l`
        filtered_rate=`echo "scale=4; 100*$filtered_reads/$raw_trimmed_reads" | bc -l`
      else
        raw_reads=`bc <<< $(grep "in total" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$raw_supplementarysecondary_reads`
        raw_trimmed_reads="NULL"
        mapped_reads_rate=`echo "scale=4; 100*$mapped_reads/$raw_reads" | bc -l`
        raw_trimmed_rate="NULL"
        filtered_rate=`echo "scale=4; 100*$filtered_reads/$raw_reads" | bc -l`
    fi
    filtered_mito_reads=$(sambamba view -F "not duplicate" -c $bam_file chrM)
    filtered_mito_rate=$(echo "scale=4; 100*$filtered_mito_reads/$filtered_mapped_reads" | bc -l)
    echo -e "$sample\\t$mark_name\\t$raw_reads\\t$raw_trimmed_reads\\t$raw_trimmed_rate\\t$mapped_reads\\t$mapped_reads_rate\\t$filtered_reads\\t$filtered_rate\\t$filtered_mapped_reads\\t$filtered_mapped_rate\\t$filtered_dup_reads\\t$filtered_dup_rate\\t$filtered_dedup_reads\\t$filtered_mito_reads\\t$filtered_mito_rate" >> {metrics_file}
  done
done && \\
sed -i -e "1 i\\Sample\\tMark Name\\tRaw Reads #\\tRemaining Reads after Trimming #\\tRemaining Reads after Trimming %\\tAligned Trimmed Reads #\\tAligned Trimmed Reads %\\tRemaining Reads after Filtering #\\tRemaining Reads after Filtering %\\tAligned Filtered Reads #\\tAligned Filtered Reads %\\tDuplicate Reads #\\tDuplicate Reads %\\tFinal Aligned Reads # without Duplicates\\tMitochondrial Reads #\\tMitochondrial Reads %" {metrics_file} && \\
mkdir -p {report_dir} && \\
cp {metrics_file} {report_metrics_file} && \\
sample_table=`LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{$1 = $1; print $0}}}}' {report_metrics_file}` && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable sample_table="$sample_table" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}
""".format(
                    sambamba=config.param('DEFAULT', 'module_sambamba'),
                    metrics_dir=metrics_output_directory,
                    metrics_file=metrics_file,
                    # samples=" ".join([sample.name for sample in self.samples]),
                    samples_associative_array=" ".join(samples_associative_array),
                    alignment_dir=self.output_dirs['alignment_output_directory'],
                    report_dir=self.output_dirs['report_output_directory'],
                    trim_metrics_file=trim_metrics_file,
                    report_metrics_file=report_metrics_file,
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                name="metrics_report",  # ".".join([sample.name for sample in self.samples]),
                samples=self.samples,
                removable_files=[report_metrics_file],
                report_files=[report_file]
            )
        )
        return jobs

    def homer_make_tag_directory(self):
        """
        The Homer Tag directories, used to check for quality metrics, are computed at this step.
        """

        jobs = []
        for sample in self.samples:
            for mark_name in sample.marks:
                alignment_file = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name,
                                              sample.name + "." + mark_name + ".sorted.dup.filtered.bam")
                output_dir = os.path.join(self.output_dirs['homer_output_directory'], sample.name,
                                          sample.name + "." + mark_name)
                other_options = config.param('homer_make_tag_directory', 'other_options', required=False)
                genome = config.param('homer_make_tag_directory', 'genome', required=False) if config.param('homer_make_tag_directory', 'genome', required=False) else self.ucsc_genome

                job = homer.makeTagDir(
                    output_dir,
                    alignment_file,
                    genome,
                    restriction_site=None,
                    illuminaPE=False,
                    other_options=other_options
                )
                job.name = "homer_make_tag_directory." + sample.name + "." + mark_name
                job.removable_files = [output_dir]
                jobs.append(job)

        return jobs

    def qc_metrics(self):
        """
        Sequencing quality metrics as tag count, tag autocorrelation, sequence bias and GC bias are generated.
        """

        # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        # if self.contrasts:
        #     design_file = os.path.relpath(self.args.design.name, self.output_dir)

        readset_file = os.path.relpath(self.args.readsets.name, self.output_dir)

        report_file = os.path.join(self.output_dirs['report_output_directory'], "ChipSeq.qc_metrics.md")
        output_files = [os.path.join(self.output_dirs['graphs_output_directory'],
                                     sample.name + "." + mark_name + "_QC_Metrics.ps") for sample in self.samples
                        for mark_name in sample.marks] + [report_file]

        jobs = []

        jobs.append(
            Job(
                [os.path.join(self.output_dirs['homer_output_directory'], sample.name,
                              sample.name + "." + mark_name, "tagInfo.txt") for sample in self.samples for mark_name
                 in sample.marks],
                output_files,
                [
                    ['qc_plots_R', 'module_mugqic_tools'],
                    ['qc_plots_R', 'module_R']
                ],
                command="""\
mkdir -p {graphs_dir} && \\
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \\
  {readset_file} \\
  {output_dir} && \\
cp {report_template_dir}/{basename_report_file} {report_file} && \\
declare -A samples_associative_array=({samples_associative_array}) && \\
for sample in ${{!samples_associative_array[@]}}
do
  for mark_name in ${{samples_associative_array[$sample]}}
  do
    cp --parents {graphs_dir}/${{sample}}.${{mark_name}}_QC_Metrics.ps {report_dir}/
    convert -rotate 90 {graphs_dir}/${{sample}}.${{mark_name}}_QC_Metrics.ps {report_dir}/graphs/${{sample}}.${{mark_name}}_QC_Metrics.png
    echo -e "----\\n\\n![QC Metrics for Sample $sample and Mark $mark_name ([download high-res image]({graphs_dir}/${{sample}}.${{mark_name}}_QC_Metrics.ps))]({graphs_dir}/${{sample}}.${{mark_name}}_QC_Metrics.png)\\n" >> {report_file}
  done
done""".format(
                    samples_associative_array=" ".join(
                        ["[\"" + sample.name + "\"]=\"" + " ".join(sample.marks.keys()) + "\"" for sample in
                         self.samples]),
                    readset_file=readset_file,
                    output_dir=self.output_dir,
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file,
                    report_dir=self.output_dirs['report_output_directory'],
                    graphs_dir=self.output_dirs['graphs_output_directory']
                ),
                name="qc_plots_R",  # ".".join([sample.name for sample in self.samples]),
                samples=self.samples,
                removable_files=output_files,
                report_files=[report_file]
            )
        )

        return jobs

    def homer_make_ucsc_file(self):
        """
        Wiggle Track Format files are generated from the aligned reads using Homer.
        The resulting files can be loaded in browsers like IGV or UCSC.
        """

        jobs = []

        for sample in self.samples:
            for mark_name in sample.marks:
                tag_dir = os.path.join(self.output_dirs['homer_output_directory'], sample.name,
                                       sample.name + "." + mark_name)
                bedgraph_dir = os.path.join(self.output_dirs['tracks_output_directory'], sample.name, mark_name)
                bedgraph_file = os.path.join(bedgraph_dir, sample.name + "." + mark_name + ".ucsc.bedGraph")
                big_wig_output = os.path.join(bedgraph_dir, "bigWig", sample.name + "." + mark_name + ".bw")

                jobs.append(
                    concat_jobs([
                        bash.mkdir(bedgraph_dir),
                        homer.makeUCSCfile(
                            tag_dir,
                            bedgraph_file
                        )
                    ],
                        name="homer_make_ucsc_file." + sample.name + "." + mark_name,
                        removable_files=[bedgraph_dir]
                    )
                )

                jobs.append(
                    concat_jobs([
                        bash.mkdir(os.path.join(bedgraph_dir, "bigWig")),
                        ucsc.bedGraphToBigWig(
                            bedgraph_file,
                            big_wig_output,
                            header=True,
                            ini_section="homer_make_ucsc_file")
                    ],
                        name="homer_make_ucsc_file_bigWig." + sample.name + "." + mark_name)
                )

        report_file = os.path.join(self.output_dirs['report_output_directory'], "ChipSeq.homer_make_ucsc_file.md")
        jobs.append(
            Job(
                [os.path.join(self.output_dirs['tracks_output_directory'], sample.name, mark_name,
                              sample.name + "." + mark_name + ".ucsc.bedGraph.gz") for sample in self.samples for
                 mark_name in sample.marks],
                [report_file],
                command="""\
mkdir -p {report_dir} && \\
zip -r {report_dir}/tracks.zip tracks/*/*/*.ucsc.bedGraph.gz && \\
cp {report_template_dir}/{basename_report_file} {report_dir}/""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file,
                    report_dir=self.output_dirs['report_output_directory']
                ),
                report_files=[report_file],
                name="homer_make_ucsc_file_report"  # ".".join([sample.name for sample in self.samples])
            )
        )

        return jobs

    def macs2_callpeak(self):
        """
        Peaks are called using the MACS2 software. Different calling strategies are used for narrow and broad peaks.
        The mfold parameter used in the model building step is estimated from a peak enrichment diagnosis run.
        The estimated mfold lower bound is 10 and the estimated upper bound can vary between 15 and 100.
        The default mfold parameter of MACS2 is [10,30].
        """

        jobs = []

        samples_associative_array = []

        for sample in self.samples:
            # samples_associative_array.append("[\"" + sample.name + "\"]=\"" + " ".join(sample.marks.keys()) + "\"")
            mark_list = []
            # if no Input file
            input_file = []
            input_file_list = [os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name,
                                            sample.name + "." + mark_name + ".sorted.dup.filtered.bam") for
                               mark_name, mark_type in sample.marks.items() if mark_type == "I"]
            if len(input_file_list) > 0:
                if len(input_file_list) > 1:
                    raise Exception("Error: Sample \"" + sample.name + "\" has more than 1 Input!")
                input_file = [input_file_list[0]]
            for mark_name, mark_type in sample.marks.items():
                if mark_type != "I":
                    mark_list.append(mark_name)

                    mark_file = [
                        os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name,
                                     sample.name + "." + mark_name + ".sorted.dup.filtered.bam")]
                    output_dir = os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name)

                    ## set macs2 variables:

                    options = "--format " + ("BAMPE" if self.run_type == "PAIRED_END" else "BAM")
                    genome_size = config.param('macs2_callpeak', 'genome_size', required=False) if config.param('macs2_callpeak', 'genome_size', required=False) else self.mappable_genome_size()
                    output_prefix_name = os.path.join(output_dir, sample.name + "." + mark_name)

                    if mark_type == "B":  # Broad region
                        options += " --broad --nomodel"
                    else:  # Narrow region
                        if input_file:
                            options += " --nomodel"
                        else:
                            options += " --fix-bimodal"

                    options += " --shift " + config.param('macs2_callpeak', 'shift') if config.param(
                        'macs2_callpeak', 'shift') else ""
                    options += " --extsize " + config.param('macs2_callpeak', 'extsize') if config.param(
                        'macs2_callpeak', 'extsize') else ""
                    options += " -p " + config.param('macs2_callpeak', 'pvalue') if config.param(
                        'macs2_callpeak', 'pvalue') else ""
                    output = []
                    output.append(os.path.join(output_dir,
                                          sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                              mark_type] + "Peak"))
                    output.append(os.path.join(output_dir,
                                 sample.name + "." + mark_name + "_peaks.xls"))


                    jobs.append(
                        concat_jobs([
                            bash.mkdir(output_dir),
                            macs2.callpeak(
                                options,
                                genome_size,
                                mark_file,
                                input_file,
                                output_prefix_name,
                                output
                            )
                        ],
                            name="macs2_callpeak." + sample.name + "." + mark_name,
                            removable_files=[output_dir]
                        )
                    )

                    ## For ihec: exchange peak score by log10 q-value and generate bigBed
                    jobs.append(
                        concat_jobs([
                            Job([os.path.join(output_dir,
                                              sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                                  mark_type] + "Peak")],
                                [os.path.join(output_dir,
                                              sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                                  mark_type] + "Peak.bed")],
                                command="""\
awk '{{if ($9 > 1000) {{$9 = 1000}}; printf( \"%s\\t%s\\t%s\\t%s\\t%0.f\\n\", $1,$2,$3,$4,$9)}}' {peak_file} > {peak_bed_file}""".format(
                                    peak_file=os.path.join(output_dir, sample.name + "." + mark_name + "_peaks." +
                                                           self.mark_type_conversion[mark_type] + "Peak"),
                                    peak_bed_file=os.path.join(output_dir,
                                                               sample.name + "." + mark_name + "_peaks." +
                                                               self.mark_type_conversion[mark_type] + "Peak.bed")
                                )
                                ),
                            ucsc.bedToBigBed(
                                os.path.join(output_dir,
                                             sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                                 mark_type] + "Peak.bed"),
                                os.path.join(output_dir,
                                             sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                                 mark_type] + "Peak.bb")
                            )
                        ],
                            name="macs2_callpeak_bigBed." + sample.name + "." + mark_name
                        )
                    )
                # Else if mark type is Input
                else:
                    log.warning("Mark " + mark_name + " for Sample " + sample.name + " is an Input ... skipping")
            samples_associative_array.append("[\"" + sample.name + "\"]=\"" + " ".join(mark_list) + "\"")

        report_file = os.path.join(self.output_dirs['report_output_directory'], "ChipSeq.macs2_callpeak.md")
        jobs.append(
            Job(
                [os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name,
                              sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                  mark_type] + "Peak") for sample in self.samples for mark_name, mark_type in
                 sample.marks.items() if mark_type != "I"],
                [report_file],
                command="""\
mkdir -p {report_dir} && \\
cp {report_template_dir}/{basename_report_file} {report_dir}/ && \\
declare -A samples_associative_array=({samples_associative_array}) && \\
for sample in ${{!samples_associative_array[@]}}
do
  for mark_name in ${{samples_associative_array[$sample]}}
  do
    cp -a --parents {macs_dir}/$sample/$mark_name/ {report_dir}/ && \\
    echo -e "* [Peak Calls File for Sample $sample and Mark $mark_name]({macs_dir}/$sample/$mark_name/${{sample}}.${{mark_name}}_peaks.xls)" >> {report_file}
  done
done""".format(
                    samples_associative_array=" ".join(samples_associative_array),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file,
                    macs_dir=self.output_dirs['macs_output_directory'],
                    report_dir=self.output_dirs['report_output_directory']
                ),
                report_files=[report_file],
                name="macs2_callpeak_report"  # ".".join([sample.name for sample in self.samples])
            )
        )

        return jobs

    def macs2_atacseq_callpeak(self):
        """
        Peaks are called using the MACS2 software. Different calling strategies are used for narrow and broad peaks.
        The mfold parameter used in the model building step is estimated from a peak enrichment diagnosis run.
        The estimated mfold lower bound is 10 and the estimated upper bound can vary between 15 and 100.
        The default mfold parameter of MACS2 is [10,30].
        """

        jobs = []

        samples_associative_array = []

        for sample in self.samples:
            # samples_associative_array.append("[\"" + sample.name + "\"]=\"" + " ".join(sample.marks.keys()) + "\"")
            mark_list = []
            # if no Input file
            input_file = []
            input_file_list = [os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name,
                                            sample.name + "." + mark_name + ".sorted.dup.filtered.bam") for
                               mark_name, mark_type in sample.marks.items() if mark_type == "I"]
            if len(input_file_list) > 0:
                if len(input_file_list) > 1:
                    raise Exception("Error: Sample \"" + sample.name + "\" has more than 1 Input!")
                input_file = [input_file_list[0]]
            for mark_name, mark_type in sample.marks.items():
                if mark_type != "I":
                    mark_list.append(mark_name)

                    mark_file = [
                        os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name,
                                     sample.name + "." + mark_name + ".sorted.dup.filtered.bam")]
                    # control_files = [os.path.join(self.output_dirs['alignment_output_directory'], sample.name, sample.name + ".sorted.dup.filtered.bam") for sample in contrast.controls]
                    output_dir = os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name)

                    ## set macs2 variables:

                    options = "--format " + ("BAMPE" if self.run_type == "PAIRED_END" else "BAM")
                    genome_size = config.param('macs2_callpeak', 'genome_size', required=False) if config.param('macs2_callpeak', 'genome_size', required=False) else self.mappable_genome_size()
                    output_prefix_name = os.path.join(output_dir, sample.name + "." + mark_name)
                    # output = os.path.join(output_dir,
                    #                       sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                    #                           mark_type] + "Peak")

                    output = []
                    output.append(os.path.join(output_dir,
                                          sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                              mark_type] + "Peak"))
                    output.append(os.path.join(output_dir,
                                 sample.name + "." + mark_name + "_peaks.xls"))
                    # other_options = " --broad --nomodel --bdg --SPMR --keep-dup all"
                    options += " --nomodel --call-summits"
                    options += " --shift " + config.param('macs2_callpeak', 'shift') if config.param(
                        'macs2_callpeak', 'shift') else " --shift -75 "
                    options += " --extsize " + config.param('macs2_callpeak', 'extsize') if config.param(
                        'macs2_callpeak', 'extsize') else " --extsize 150 "
                    options += " -p " + config.param('macs2_callpeak', 'pvalue') if config.param(
                        'macs2_callpeak', 'pvalue') else " -p 0.01 "

                    jobs.append(
                        concat_jobs([
                            bash.mkdir(output_dir),
                            macs2.callpeak(
                                options,
                                genome_size,
                                mark_file,
                                input_file,
                                output_prefix_name,
                                output
                            )
                        ],
                            name="macs2_callpeak." + sample.name + "." + mark_name,
                            removable_files=[output_dir]
                        )
                    )

                    ## For ihec: exchange peak score by log10 q-value and generate bigBed
                    jobs.append(
                        concat_jobs([
                            Job([os.path.join(output_dir,
                                              sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                                  mark_type] + "Peak")],
                                [os.path.join(output_dir,
                                              sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                                  mark_type] + "Peak.bed")],
                                command="""\
awk '{{if ($9 > 1000) {{$9 = 1000}}; printf( \"%s\\t%s\\t%s\\t%s\\t%0.f\\n\", $1,$2,$3,$4,$9)}}' {peak_file} > {peak_bed_file}""".format(
                                    peak_file=os.path.join(output_dir, sample.name + "." + mark_name + "_peaks." +
                                                           self.mark_type_conversion[mark_type] + "Peak"),
                                    peak_bed_file=os.path.join(output_dir,
                                                               sample.name + "." + mark_name + "_peaks." +
                                                               self.mark_type_conversion[mark_type] + "Peak.bed")
                                )
                                ),
                            ucsc.bedToBigBed(
                                os.path.join(output_dir,
                                             sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                                 mark_type] + "Peak.bed"),
                                os.path.join(output_dir,
                                             sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                                 mark_type] + "Peak.bb")
                            )
                        ],
                            name="macs2_callpeak_bigBed." + sample.name + "." + mark_name
                        )
                    )
                # Else if mark type is Input
                else:
                    log.warning("Mark " + mark_name + " for Sample " + sample.name + " is an Input ... skipping")
            samples_associative_array.append("[\"" + sample.name + "\"]=\"" + " ".join(mark_list) + "\"")

        report_file = os.path.join(self.output_dirs['report_output_directory'], "ChipSeq.macs2_callpeak.md")
        jobs.append(
            Job(
                [os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name,
                              sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                  mark_type] + "Peak") for sample in self.samples for mark_name, mark_type in
                 sample.marks.items() if mark_type != "I"],
                [report_file],
                command="""\
mkdir -p {report_dir} && \\
cp {report_template_dir}/{basename_report_file} {report_dir}/ && \\
declare -A samples_associative_array=({samples_associative_array}) && \\
for sample in ${{!samples_associative_array[@]}}
do
  for mark_name in ${{samples_associative_array[$sample]}}
  do
    cp -a --parents {macs_dir}/$sample/$mark_name/ {report_dir}/ && \\
    echo -e "* [Peak Calls File for Sample $sample and Mark $mark_name]({macs_dir}/$sample/$mark_name/${{sample}}.${{mark_name}}_peaks.xls)" >> {report_file}
  done
done""".format(
                    samples_associative_array=" ".join(samples_associative_array),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file,
                    macs_dir=self.output_dirs['macs_output_directory'],
                    report_dir=self.output_dirs['report_output_directory']
                ),
                report_files=[report_file],
                name="macs2_callpeak_report"  # ".".join([sample.name for sample in self.samples])
            )
        )
        return jobs

    def differential_binding(self):
        """
        Performs differential binding analysis using [DiffBind](http://bioconductor.org/packages/release/bioc/html/DESeq.html)
        Merge the results of the analysis in a single csv file.
        html report will be generated to QC samples and check how well differential binding analysis was performed.
        """
        jobs = []
        minOverlap = config.param('differential_binding', 'minOverlap')
        minMembers = config.param('differential_binding', 'minMembers')
        method = config.param('differential_binding', 'method')
        # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        readset_file = os.path.relpath(self.args.readsets.name, self.output_dir)
        if self.contrasts:
            design_file = os.path.relpath(self.args.design.name, self.output_dir)
        else:
            log.info("Comparison column is not defined. Skipping differential binding analysis...")
        mark_list = []

        #if control samples and treatment samples are less than one diff analysis will not be executed
        for contrast in self.contrasts:
            bam_list = []
            controls_count = len(contrast.controls)
            treatments_count = len(contrast.treatments)
            if controls_count < 2 or treatments_count < 2:
                log.info(
                    "At leaset two treatments and  controls should be defined. Skipping differential binding analysis for "+contrast.name +" ...")
            else:
                for control in contrast.controls:
                    control_sample_name, control_mark_name = control.split("-.-")
                    for sample in self.samples:
                        input_file_list = [
                            os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name,
                                         sample.name + "." + mark_name + ".sorted.dup.filtered.bam") for mark_name, mark_type in sample.marks.items() if
                            mark_type == "I" and sample.name == control_sample_name]
                        bam_list.append(input_file_list)

                        input_file_list = [
                            os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name,
                                         sample.name + "." + mark_name + ".sorted.dup.filtered.bam") for
                            mark_name, mark_type in sample.marks.items() if
                            mark_type != "I" and sample.name == control_sample_name and mark_name ==
                            control_mark_name]
                        bam_list.append(input_file_list)

                        input_file_list = [
                            os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name,
                                        sample.name + "." + mark_name + "_peaks.xls") for
                            # os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name, sample.name + "." + mark_name + "_peaks." +
                            #                          self.mark_type_conversion[mark_type] + "Peak") for
                            mark_name, mark_type in sample.marks.items() if
                            mark_type != "I" and sample.name == control_sample_name and mark_name ==
                            control_mark_name]

                        bam_list.append(input_file_list)

                for treatment in contrast.treatments:
                    treatment_sample_name, treatment_mark_name = treatment.split("-.-")
                    for sample in self.samples:
                        input_file_list = [
                            os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name,
                                         sample.name + "." + mark_name + ".sorted.dup.filtered.bam") for
                            mark_name, mark_type in sample.marks.items() if
                            mark_type == "I" and sample.name == treatment_sample_name]
                        bam_list.append(input_file_list)

                        input_file_list = [
                            os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name,
                                         sample.name + "." + mark_name + ".sorted.dup.filtered.bam") for
                            mark_name, mark_type in sample.marks.items() if
                            mark_type != "I" and sample.name == treatment_sample_name and mark_name ==
                            treatment_mark_name]
                        bam_list.append(input_file_list)

                        input_file_list = [
                            os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name,
                                        sample.name + "." + mark_name + "_peaks.xls") for
                            #os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name, sample.name + "." + mark_name + "_peaks." +
                            #                         self.mark_type_conversion[mark_type] + "Peak") for
                            mark_name, mark_type in sample.marks.items() if
                            mark_type != "I" and sample.name == treatment_sample_name and mark_name ==
                            treatment_mark_name]

                        bam_list.append(input_file_list)
                bam_list = filter(None, bam_list)
                bam_list = [item for sublist in bam_list for item in sublist]
                diffbind_job = differential_binding.diffbind(bam_list, contrast.name, design_file, readset_file,
                                                             self.output_dirs['dba_output_directory'],
                                                             self.output_dirs['alignment_output_directory'],
                                                             self.output_dirs['macs_output_directory'], minOverlap,
                                                             minMembers, method)
                diffbind_job.samples = self.samples
                diffbind_job.name = "_".join(("differential_binding.diffbind.contrast", contrast.name))
                jobs.append(diffbind_job)

        return jobs

    def homer_annotate_peaks(self):
            """
            The peaks called previously are annotated with HOMER using RefSeq annotations for the reference genome.
            Gene ontology and genome ontology analysis are also performed at this stage.
            """

            jobs = []

            samples_associative_array = []

            for sample in self.samples:
                # samples_associative_array.append("[\"" + sample.name + "\"]=\"" + " ".join(sample.marks.keys()) + "\"")
                mark_list = []
                for mark_name, mark_type in sample.marks.items():
                    if mark_type != "I":
                        mark_list.append(mark_name)

                        peak_file = os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name,
                                                 sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                                     mark_type] + "Peak")
                        output_dir = os.path.join(self.output_dirs['anno_output_directory'], sample.name, mark_name)
                        output_prefix = os.path.join(output_dir, sample.name + "." + mark_name)
                        annotation_file = output_prefix + ".annotated.csv"

                        genome = config.param('homer_annotate_peaks', 'genome', required=False) if config.param('homer_annotate_peaks', 'genome', required=False) else self.ucsc_genome
                        genome_size = config.param('homer_annotate_peaks', 'genome_size', required=False) if config.param('homer_annotate_peaks', 'genome_size', required=False) else self.mappable_genome_size()

                        jobs.append(
                            concat_jobs([
                                bash.mkdir(output_dir),
                                homer.annotatePeaks(
                                    peak_file,
                                    genome,
                                    output_dir,
                                    annotation_file,
                                    genome_size
                                ),
                                Job(
                                    [annotation_file],
                                    [
                                        output_prefix + ".tss.stats.csv",
                                        output_prefix + ".exon.stats.csv",
                                        output_prefix + ".intron.stats.csv",
                                        output_prefix + ".tss.distance.csv"
                                    ],
                                    [['homer_annotate_peaks', 'module_perl'],
                                     ['homer_annotate_peaks', 'module_mugqic_tools']],
                                    command="""\
    perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
      "{annotation_file}",
      "{output_prefix}",
      {proximal_distance},
      {distal_distance},
      {distance5d_lower},
      {distance5d_upper},
      {gene_desert_size}
    )'""".format(
                                        annotation_file=annotation_file,
                                        output_prefix=output_prefix,
                                        proximal_distance=config.param('homer_annotate_peaks', 'proximal_distance',
                                                                       param_type='int'),
                                        distal_distance=config.param('homer_annotate_peaks', 'distal_distance',
                                                                     param_type='int'),
                                        distance5d_lower=config.param('homer_annotate_peaks', 'distance5d_lower',
                                                                      param_type='int'),
                                        distance5d_upper=config.param('homer_annotate_peaks', 'distance5d_upper',
                                                                      param_type='int'),
                                        gene_desert_size=config.param('homer_annotate_peaks', 'gene_desert_size',
                                                                      param_type='int')
                                    ),
                                    removable_files=[
                                        os.path.join(self.output_dirs['anno_output_directory'], sample.name,
                                                     mark_name)],
                                )
                            ],
                                name="homer_annotate_peaks." + sample.name + "." + mark_name)
                        )

                    else:
                        log.warning("Mark " + mark_name + " for Sample " + sample.name + " is an Input ... skipping")
                samples_associative_array.append("[\"" + sample.name + "\"]=\"" + " ".join(mark_list) + "\"")

            report_file = os.path.join(self.output_dirs['report_output_directory'], "ChipSeq.homer_annotate_peaks.md")
            jobs.append(
                Job(
                    [os.path.join(self.output_dirs['anno_output_directory'], sample.name, mark_name,
                                  sample.name + "." + mark_name + ".annotated.csv") for sample in self.samples for
                     mark_name, mark_type in sample.marks.items() if mark_type != "I"],
                    [report_file],
                    command="""\
    mkdir -p {report_dir}/annotation/ && \\
    cp {report_template_dir}/{basename_report_file} {report_dir} && \\
    declare -A samples_associative_array=({samples_associative_array}) && \\
    for sample in ${{!samples_associative_array[@]}}
    do
      for mark_name in ${{samples_associative_array[$sample]}}
      do
        rsync -rvP annotation/$sample {report_dir}/annotation/ && \\
        echo -e "* [Gene Annotations for Sample $sample and Mark $mark_name](annotation/$sample/$mark_name/${{sample}}.${{mark_name}}.annotated.csv)\\n* [HOMER Gene Ontology Annotations for Sample $sample and Mark $mark_name](annotation/$sample/$mark_name/geneOntology.html)\\n* [HOMER Genome Ontology Annotations for Sample $sample and Mark $mark_name](annotation/$sample/$mark_name/GenomeOntology.html)" >> {report_file}
      done
    done""".format(
                        samples_associative_array=" ".join(samples_associative_array),
                        report_template_dir=self.report_template_dir,
                        basename_report_file=os.path.basename(report_file),
                        report_file=report_file,
                        report_dir=self.output_dirs['report_output_directory']
                    ),
                    report_files=[report_file],
                    name="homer_annotate_peaks_report"  # ".".join([sample.name for sample in self.samples])
                )
            )

            return jobs

    def homer_find_motifs_genome(self):
            """
            De novo and known motif analysis per design are performed using HOMER.
            """

            jobs = []

            counter = 0

            samples_associative_array = []

            for sample in self.samples:
                mark_list = []
                # samples_associative_array.append("[\"" + sample.name + "\"]=\"" + " ".join(sample.marks.keys()) + "\"")
                for mark_name, mark_type in sample.marks.items():
                    # Don't find motifs for broad peaks
                    if mark_type == "N":
                        mark_list.append(mark_name)

                        peak_file = os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name,
                                                 sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                                     mark_type] + "Peak")
                        output_dir = os.path.join(self.output_dirs['anno_output_directory'], sample.name, mark_name)

                        genome = config.param('homer_annotate_peaks', 'genome', required=False) if config.param('homer_annotate_peaks', 'genome', required=False) else self.ucsc_genome

                        jobs.append(
                            concat_jobs([
                                bash.mkdir(output_dir),
                                homer.findMotifsGenome(
                                    peak_file,
                                    genome,
                                    output_dir,
                                    config.param('homer_find_motifs_genome', 'threads', param_type='posint')
                                )
                            ],
                                name="homer_find_motifs_genome." + sample.name + "." + mark_name,
                                removable_files=[
                                    os.path.join(self.output_dirs['anno_output_directory'], sample.name, mark_name)]
                            )
                        )
                        counter = counter + 1
                    else:
                        # log.warning("No treatment found for contrast " + contrast.name + "... skipping")
                        # log.warning("Contrast " + contrast.name + " is broad; homer_find_motifs_genome is run on narrow peaks ... skipping")
                        log.warning(
                            "Mark " + mark_name + " for Sample " + sample.name + " is not Narrow; homer_find_motifs_genome is run on narrow peaks ... skipping")
                samples_associative_array.append("[\"" + sample.name + "\"]=\"" + " ".join(mark_list) + "\"")

            if counter > 0:
                report_file = os.path.join(self.output_dirs['report_output_directory'],
                                           "ChipSeq.homer_find_motifs_genome.md")
                jobs.append(
                    Job(
                        [os.path.join(self.output_dirs['anno_output_directory'], sample.name, mark_name,
                                      "homerResults.html") for sample in self.samples for mark_name, mark_type in
                         sample.marks.items() if mark_type == "N"] +
                        [os.path.join(self.output_dirs['anno_output_directory'], sample.name, mark_name,
                                      "knownResults.html") for sample in self.samples for mark_name, mark_type in
                         sample.marks.items() if mark_type == "N"],
                        [report_file],
                        command="""\
    mkdir -p {report_dir}/annotation/ && \\
    cp {report_template_dir}/{basename_report_file} {report_dir}/ && \\
    declare -A samples_associative_array=({samples_associative_array}) && \\
    for sample in ${{!samples_associative_array[@]}}
    do
      for mark_name in ${{samples_associative_array[$sample]}}
      do
        rsync -rvP annotation/$sample {report_dir}/annotation/ && \\
        echo -e "* [HOMER _De Novo_ Motif Results for Sample $sample and Mark $mark_name](annotation/$sample/$mark_name/homerResults.html)\\n* [HOMER Known Motif Results for Sample $sample and Mark $mark_name](annotation/$sample/$mark_name/knownResults.html)" >> {report_file}
      done
    done""".format(
                            samples_associative_array=" ".join(samples_associative_array),
                            report_template_dir=self.report_template_dir,
                            basename_report_file=os.path.basename(report_file),
                            report_file=report_file,
                            report_dir=self.output_dirs['report_output_directory']
                        ),
                        report_files=[report_file],
                        name="homer_find_motifs_genome_report"  # ".".join([sample.name for sample in self.samples])
                    )
                )

            return jobs

    def annotation_graphs(self):
            """
            The peak location statistics. The following peak location statistics are generated per design:
            proportions of the genomic locations of the peaks. The locations are: Gene (exon or intron),
            Proximal ([0;2] kb upstream of a transcription start site), Distal ([2;10] kb upstream
            of a transcription start site), 5d ([10;100] kb upstream of a transcription start site),
            Gene desert (>= 100 kb upstream or downstream of a transcription start site), Other (anything
            not included in the above categories); The distribution of peaks found within exons and introns;
            The distribution of peak distance relative to the transcription start sites (TSS);
            the Location of peaks per design.
            """

            # If --design <design_file> option is missing, self.contrasts call will raise an Exception
            # if self.contrasts:
            #     design_file = os.path.relpath(self.args.design.name, self.output_dir)

            readset_file = os.path.relpath(self.args.readsets.name, self.output_dir)

            input_files = []
            output_files = []
            samples_associative_array = []
            for sample in self.samples:
                mark_list = []
                for mark_name, mark_type in sample.marks.items():
                    if mark_type == "N":
                        annotation_prefix = os.path.join(self.output_dirs['anno_output_directory'], sample.name,
                                                         mark_name, sample.name + "." + mark_name)
                        input_files.append(annotation_prefix + ".tss.stats.csv")
                        input_files.append(annotation_prefix + ".exon.stats.csv")
                        input_files.append(annotation_prefix + ".intron.stats.csv")
                        input_files.append(annotation_prefix + ".tss.distance.csv")
                        mark_list.append(mark_name)
                samples_associative_array.append("[\"" + sample.name + "\"]=\"" + " ".join(mark_list) + "\"")

                peak_stats_file = os.path.join(self.output_dirs['anno_output_directory'], sample.name, "peak_stats.csv")
                output_files.append(peak_stats_file)
            report_file = os.path.join(self.output_dirs['report_output_directory'], "ChipSeq.annotation_graphs.md")
            output_files.append(report_file)

            jobs = []

            jobs.append(
                Job(
                    input_files,
                    output_files,
                    [
                        ['annotation_graphs', 'module_mugqic_tools'],
                        ['annotation_graphs', 'module_R'],
                        ['annotation_graphs', 'module_pandoc']
                    ],
                    command="""\
    cp /dev/null annotation/peak_stats_AllSamples.csv && \\
    mkdir -p {graphs_dir} && \\
    Rscript $R_TOOLS/chipSeqgenerateAnnotationGraphs.R \\
      {readset_file} \\
      {output_dir} && \\
    declare -A samples_associative_array=({samples_associative_array}) && \\
    for sample in ${{!samples_associative_array[@]}}
    do
        header=$(head -n 1 annotation/$sample/peak_stats.csv)
        tail -n+2 annotation/$sample/peak_stats.csv >> annotation/peak_stats_AllSamples.csv
    done && \\
    sed -i -e "1 i\\\$header" annotation/peak_stats_AllSamples.csv && \\
    mkdir -p {report_dir}/annotation/$sample && \\
    cp annotation/peak_stats_AllSamples.csv {report_dir}/annotation/peak_stats_AllSamples.csv && \\
    peak_stats_table=`LC_NUMERIC=en_CA awk -F "," '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{print $1, $2,  sprintf("%\\47d", $3), $4, sprintf("%\\47.1f", $5), sprintf("%\\47.1f", $6), sprintf("%\\47.1f", $7), sprintf("%\\47.1f", $8)}}}}' annotation/peak_stats_AllSamples.csv`
    pandoc --to=markdown \\
        --template {report_template_dir}/{basename_report_file} \\
        --variable peak_stats_table="$peak_stats_table" \\
        --variable proximal_distance="{proximal_distance}" \\
        --variable distal_distance="{distal_distance}" \\
        --variable distance5d_lower="{distance5d_lower}" \\
        --variable distance5d_upper="{distance5d_upper}" \\
        --variable gene_desert_size="{gene_desert_size}" \\
        {report_template_dir}/{basename_report_file} \\
        > {report_file} && \\
    for sample in ${{!samples_associative_array[@]}}
    do
      cp annotation/$sample/peak_stats.csv {report_dir}/annotation/$sample/peak_stats.csv && \\
      for mark_name in ${{samples_associative_array[$sample]}}
      do
        cp --parents {graphs_dir}/${{sample}}.${{mark_name}}_Misc_Graphs.ps {report_dir}/
        convert -rotate 90 {graphs_dir}/${{sample}}.${{mark_name}}_Misc_Graphs.ps {report_dir}/graphs/${{sample}}.${{mark_name}}_Misc_Graphs.png
        echo -e "----\\n\\n![Annotation Statistics for Sample $sample and Mark $mark_name ([download high-res image]({graphs_dir}/${{sample}}.${{mark_name}}_Misc_Graphs.ps))]({graphs_dir}/${{sample}}.${{mark_name}}_Misc_Graphs.png)\\n" >> {report_file}
      done
    done""".format(
                        readset_file=readset_file,
                        output_dir=self.output_dir,
                        # peak_stats_file=peak_stats_file,
                        # samples_associative_array=" ".join(["[\"" + sample.name + "\"]=\"" + " ".join(sample.marks.keys()) + "\"" for sample in self.samples]),
                        samples_associative_array=" ".join(samples_associative_array),
                        # contrasts=" ".join([contrast.real_name for contrast in self.contrasts if contrast.type == 'narrow' and contrast.treatments]),
                        proximal_distance=config.param('homer_annotate_peaks', 'proximal_distance', param_type='int') / -1000,
                        distal_distance=config.param('homer_annotate_peaks', 'distal_distance', param_type='int') / -1000,
                        distance5d_lower=config.param('homer_annotate_peaks', 'distance5d_lower', param_type='int') / -1000,
                        distance5d_upper=config.param('homer_annotate_peaks', 'distance5d_upper', param_type='int') / -1000,
                        gene_desert_size=config.param('homer_annotate_peaks', 'gene_desert_size', param_type='int') / 1000,
                        report_template_dir=self.report_template_dir,
                        basename_report_file=os.path.basename(report_file),
                        report_file=report_file,
                        report_dir=self.output_dirs['report_output_directory'],
                        graphs_dir=self.output_dirs['graphs_output_directory'],
                        merged_peak_stats="peak_stats_AllSamples.csv"
                    ),
                    name="annotation_graphs",  # ".".join([sample.name for sample in self.samples]),
                    report_files=[report_file],
                    removable_files=output_files
                )
            )

            return jobs


    def run_spp(self):
            """
            runs spp to estimate NSC and RSC ENCODE metrics. For more information: https://github.com/kundajelab/phantompeakqualtools
            """
            jobs = []

            for sample in self.samples:
                for mark_name in sample.marks:
                    alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], sample.name,
                                                       mark_name)
                    sample_merge_mdup_bam = os.path.join(alignment_directory,
                                                         sample.name + "." + mark_name + ".sorted.dup.filtered.bam")
                    output_dir = os.path.join(self.output_dirs['ihecM_output_directory'], sample.name, mark_name)
                    output = os.path.join(output_dir, sample.name + "." + mark_name + ".crosscor")

                    jobs.append(
                        concat_jobs([
                            bash.mkdir(output_dir),
                            Job(
                                [sample_merge_mdup_bam],
                                [output],
                                [
                                    ['run_spp', 'module_samtools'],
                                    ['run_spp', 'module_mugqic_tools'],
                                    ['run_spp', 'module_R']
                                ],
                                command="""\
    cat /dev/null > {output} && \\
    Rscript $R_TOOLS/run_spp.R -c={sample_merge_mdup_bam} -savp -out={output} -rf -tmpdir={tmp_dir}""".format(
                                    sample_merge_mdup_bam=sample_merge_mdup_bam,
                                    output=output,
                                    tmp_dir=config.param('run_spp', 'tmp_dir')
                                )
                            )
                        ],
                            name="run_spp." + sample.name + "." + mark_name)
                    )

            jobs.append(
                Job(
                    [os.path.join(self.output_dirs['ihecM_output_directory'], sample.name, mark_name,
                                  sample.name + "." + mark_name + ".crosscor") for sample in self.samples for
                     mark_name, mark_type in sample.marks.items()],
                    [os.path.join(self.output_dirs['ihecM_output_directory'], sample.name, sample.name + ".crosscor")
                     for sample in self.samples],
                    [],
                    command="""\
    declare -A samples_associative_array=({samples_associative_array}) && \\
    for sample in ${{!samples_associative_array[@]}}
    do
      echo -e "Filename\\tnumReads\\testFragLen\\tcorr_estFragLen\\tPhantomPeak\\tcorr_phantomPeak\\targmin_corr\\tmin_corr\\tNormalized SCC (NSC)\\tRelative SCC (RSC)\\tQualityTag)" > ihec_metrics/${{sample}}/${{sample}}.crosscor
      for mark_name in ${{samples_associative_array[$sample]}}
      do
        cat ihec_metrics/${{sample}}/${{mark_name}}/${{sample}}.${{mark_name}}.crosscor >> ihec_metrics/${{sample}}/${{sample}}.crosscor
      done
    done""".format(
                        samples_associative_array=" ".join(
                            ["[\"" + sample.name + "\"]=\"" + " ".join(sample.marks.keys()) + "\"" for sample in
                             self.samples])
                    ),
                    name="run_spp_report"  # ".".join([sample.name for sample in self.samples])
                )
            )

            return jobs

    def ihec_metrics(self):
            """
            Generate IHEC's standard metrics.
            """
            jobs = []

            alignment_dir = self.output_dirs['alignment_output_directory']
            # output_dir = self.output_dirs['ihecM_output_directory']

            # samples_associative_array = []
            metrics_to_merge = []

            for sample in self.samples:
                mark_list = []
                # if no Input file
                input_file = {}
                input_file_list = [mark_name for mark_name, mark_type in sample.marks.items() if mark_type == "I"]
                if len(input_file_list) > 0:
                    if len(input_file_list) > 1:
                        raise Exception("Error: Sample \"" + sample.name + "\" has more than 1 Input!")
                    input_file[sample.name] = input_file_list[0]
                for mark_name, mark_type in sample.marks.items():
                    if mark_type != "I":
                        mark_list.append(mark_name)

                        chip_bam = os.path.join(alignment_dir, sample.name, mark_name,
                                                sample.name + "." + mark_name + ".sorted.dup.bam")
                        chip_bed = os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name,
                                                sample.name + "." + mark_name + "_peaks." + self.mark_type_conversion[
                                                    mark_type] + "Peak.bed")
                        output_dir = os.path.join(self.output_dirs['ihecM_output_directory'], sample.name)
                        crosscor_input = os.path.join(self.output_dirs['ihecM_output_directory'], sample.name,
                                                      sample.name + ".crosscor")
                        genome = config.param('IHEC_chipseq_metrics', 'assembly')

                        if not input_file:
                            input_name = "no_input"
                            input_bam = None
                        else:
                            input_name = input_file[sample.name]  # "".join(input_file.keys())
                            input_bam = os.path.join(alignment_dir, sample.name, input_name,
                                                     sample.name + "." + input_name + ".sorted.dup.bam")  # input_file[sample.name]

                        jobs.append(
                            concat_jobs([
                                bash.mkdir(output_dir),
                                tools.sh_ihec_chip_metrics(
                                    chip_bam=chip_bam,
                                    input_bam=input_bam,
                                    sample_name=sample.name,
                                    input_name=input_name,
                                    chip_name=mark_name,
                                    chip_type=self.mark_type_conversion[mark_type],
                                    chip_bed=chip_bed,
                                    output_dir=output_dir,
                                    assembly=genome,
                                    crosscor_input=crosscor_input
                                )
                            ],
                                name="IHEC_chipseq_metrics." + sample.name + "." + mark_name,
                                removable_files=[output_dir]
                            )
                        )
                        metrics_to_merge.append(os.path.join(output_dir, mark_name,
                                                             "IHEC_chipseq_metrics." + sample.name + "." + mark_name + ".tsv"))

            metrics_merged = "IHEC_chipseq_metrics_AllSamples.tsv"
            metrics_merged_out = os.path.join(self.output_dirs['ihecM_output_directory'], metrics_merged)
            report_file = os.path.join("report", "ChipSeq.ihec_metrics.md")

            jobs.append(
                Job(
                    input_files=metrics_to_merge,
                    output_files=[metrics_merged_out],
                    name="merge_ihec_metrics",  # ".".join([sample.name for sample in self.samples]),
                    command="""\
    cp /dev/null {metrics_merged} && \\
    for sample in {samples}
    do
        header=$(head -n 1 $sample | cut -f -3,5-17,30-33,35,37,39-)
        tail -n 1 $sample | cut -f -3,5-17,30-33,35,37,39- >> {metrics_merged}
    done && \\
    sample_name=`tail -n 1 $sample | cut -f 1` && \\
    input_name=`tail -n 1 $sample | cut -f 4` && \\
    input_chip_type="NA" && \\
    genome_assembly=`tail -n 1 $sample | cut -f 5` && \\
    input_core=`tail -n 1 $sample | cut -f 18-29` && \\
    input_nsc=`tail -n 1 $sample | cut -f 34` && \\
    input_rsc=`tail -n 1 $sample | cut -f 36` && \\
    input_quality=`tail -n 1 $sample | cut -f 38` && \\
    if [[ $input_name != "no_input" ]]
      then
        echo -e "${{sample_name}}\\t${{input_name}}\\t${{input_chip_type}}\\t${{genome_assembly}}\\t${{input_core}}\\tNA\\tNA\\tNA\\t${{input_nsc}}\\t${{input_rsc}}\\t${{input_quality}}\\tNA\\tNA" >> {metrics_merged}
    fi && \\
    sed -i -e "1 i\\\$header" {metrics_merged}""".format(
                        samples=" ".join(metrics_to_merge),
                        metrics_merged=metrics_merged_out
                    ),
                )
            )

            jobs.append(
                Job(
                    input_files=[metrics_merged_out],
                    output_files=[report_file],
                    # name="merge_ihec_metrics_report." + ".".join([sample.name for sample in self.samples]),
                    name="merge_ihec_metrics_report",
                    module_entries=[['merge_ihec_metrics_report', 'module_pandoc']],
                    command="""\
    mkdir -p {report_dir} && \\
    cp {metrics_merged_out} {report_dir}/{ihec_metrics_merged_table} && \\
    pandoc --to=markdown \\
      --template {report_template_dir}/{basename_report_file} \\
      --variable ihec_metrics_merged_table="{ihec_metrics_merged_table}" \\
      {report_template_dir}/{basename_report_file} \\
      > {report_file}""".format(
                        metrics_merged_out=metrics_merged_out,
                        ihec_metrics_merged_table=metrics_merged,
                        report_template_dir=self.report_template_dir,
                        basename_report_file=os.path.basename(report_file),
                        report_file=report_file,
                        report_dir=self.output_dirs['report_output_directory']
                    ),
                    report_files=[report_file]
                )
            )

            return jobs

    def multiqc_report(self):
            """
            A quality control report for all samples is generated.
            For more detailed information about the MultiQc visit: [MultiQc] (http://multiqc.info/)
            """
            ## set multiQc config file so we can customize one for every pipeline:
            jobs = []
            # yamlFile = os.path.expandvars(config.param('multiqc_report', 'MULTIQC_CONFIG_PATH'))
            input_files = []
            metrics_output_directory = self.output_dirs['metrics_output_directory']
            for sample in self.samples:
                for mark_name in sample.marks:
                    picard_prefix = os.path.join(metrics_output_directory, sample.name, mark_name,
                                                 sample.name + "." + mark_name + ".sorted.dup.filtered.all.metrics.")
                    if self.run_type == 'SINGLE_END':
                        picard_files = [
                            picard_prefix + "quality_by_cycle.pdf",
                            picard_prefix + "alignment_summary_metrics",
                            picard_prefix + "quality_by_cycle_metrics",
                            picard_prefix + "quality_distribution_metrics",
                            picard_prefix + "quality_distribution.pdf"
                        ]
                    elif self.run_type == 'PAIRED_END':
                        picard_files = [
                            picard_prefix + "base_distribution_by_cycle.pdf",
                            picard_prefix + "alignment_summary_metrics",
                            picard_prefix + "insert_size_histogram.pdf",
                            picard_prefix + "insert_size_metrics",
                            picard_prefix + "quality_by_cycle_metrics",
                            picard_prefix + "quality_by_cycle.pdf",
                            picard_prefix + "quality_distribution_metrics",
                            picard_prefix + "quality_distribution.pdf"

                        ]
                    input_files.extend(picard_files)
                    input_files.append(os.path.join(metrics_output_directory, sample.name, mark_name,
                                                    sample.name + "." + mark_name + ".sorted.dup.filtered.flagstat"))
                    homer_prefix = os.path.join(self.output_dirs['homer_output_directory'], sample.name,
                                                sample.name + "." + mark_name)
                    homer_files = [
                        os.path.join(homer_prefix, "tagGCcontent.txt"),
                        os.path.join(homer_prefix, "genomeGCcontent.txt"),
                        os.path.join(homer_prefix, "tagLengthDistribution.txt"),
                        os.path.join(homer_prefix, "tagInfo.txt")
                    ]
                    input_files.extend(homer_files)
            # input_files = [os.path.join(self.output_dirs['homer_output_directory'], sample.name, "tagInfo.txt") for sample in self.samples]
            output = os.path.join(self.output_dirs['report_output_directory'], "multiqc_report")
            log.info(output)

            job = multiqc.run(
                input_files,
                output,
                ini_section='multiqc_report'
            )
            job.name = "multiqc_report"  # ".".join([sample.name for sample in self.samples])

            jobs.append(job)

            return jobs

    def cram_output(self):
        """
        Generate long term storage version of the final alignment files in CRAM format
        Using this function will include the orginal final bam file into the  removable file list
        """

        jobs = []

        for sample in self.samples:
            for mark_name in sample.marks:
                input_bam = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name,
                                         sample.name + "." + mark_name + ".sorted.dup.filtered.bam")
                output_cram = re.sub("\.bam$", ".cram", input_bam)

                # Run samtools
                job = samtools.view(
                    input_bam,
                    output_cram,
                    options=config.param('samtools_cram_output', 'options'),
                    removable=False
                )
                job.name = "cram_output." + sample.name + "." + mark_name
                job.removable_files = input_bam

                jobs.append(job)

        return jobs

    def gatk_haplotype_caller(self):
        """
        GATK haplotype caller for snps and small indels.
        """

        jobs = []


        for sample in self.samples:
            for mark_name, mark_type in sample.marks.items():
                if not mark_type == "I":
                    alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], sample.name,
                                               mark_name)
                    haplotype_directory = os.path.join(alignment_directory, "rawHaplotypeCaller")

                    macs_output_dir = os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name)

                    input_bam = os.path.join(alignment_directory,
                                                 sample.name + "." + mark_name + ".sorted.dup.filtered.bam")

                    #peak calling bed file from MACS2 is given here to restrict the variant calling to peaks regions
                    interval_list = None
                    if mark_type == "N":
                        interval_list = os.path.join(macs_output_dir, sample.name + "." + mark_name +"_peaks.narrowPeak.bed")
                    elif mark_type == "B":
                        interval_list = os.path.join(macs_output_dir,
                                                     sample.name + "." + mark_name + "_peaks.broadPeak.bed")


                    mkdir_job = bash.mkdir(
                                haplotype_directory,
                                remove=True
                    )
                    interval_padding = config.param('gatk_haplotype_caller', 'interval_padding')
                    jobs.append(
                    concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    gatk4.haplotype_caller(
                        input_bam,
                        os.path.join(haplotype_directory, sample.name +"."+ mark_name +".hc.g.vcf.gz"),
                        interval_list=interval_list

                    )
                ],
                    name="gatk_haplotype_caller." + sample.name +"_"+ mark_name,
                    samples=[sample]
                    )
                    )

        return jobs

    def merge_and_call_individual_gvcf(self):
        """
        Merges the gvcfs of haplotype caller and also generates a per sample vcf containing genotypes.
        """

        jobs = []

        for sample in self.samples:
            for mark_name, mark_type in sample.marks.items():
                if not mark_type == "I":
                    alignment_directory = os.path.join("alignment", sample.name, mark_name)
                    haplotype_directory = os.path.join(alignment_directory, "rawHaplotypeCaller")
                    haplotype_file_prefix = os.path.join(haplotype_directory , sample.name +"." +mark_name)
                    output_haplotype_file_prefix = os.path.join("alignment", sample.name, mark_name, sample.name +"." + mark_name)

                    jobs.append(
                        concat_jobs([
                            bash.ln(
                                haplotype_file_prefix + ".hc.g.vcf.gz",
                                output_haplotype_file_prefix + ".hc.g.vcf.gz",
                                self.output_dir
                            ),
                            bash.ln(
                                haplotype_file_prefix + ".hc.g.vcf.gz.tbi",
                                output_haplotype_file_prefix + ".hc.g.vcf.gz.tbi",
                                self.output_dir
                            ),
                            gatk4.genotype_gvcf(
                                output_haplotype_file_prefix + ".hc.g.vcf.gz",
                                output_haplotype_file_prefix + ".hc.vcf.gz",
                                config.param('gatk_genotype_gvcf', 'options')
                            )
                        ],
                            name="merge_and_call_individual_gvcf.call." + sample.name +"_"+ mark_name,
                            samples=[sample]
                        )
                    )
        return jobs

    @property
    def steps(self):
        return [
            [
                self.picard_sam_to_fastq,
                self.trimmomatic,
                self.merge_trimmomatic_stats,
                self.mapping_bwa_mem_sambamba,
                self.sambamba_merge_bam_files, #5
                self.sambamba_mark_duplicates,
                self.sambamba_view_filter,
                # self.picard_mark_duplicates,
                self.metrics,
                self.homer_make_tag_directory,
                self.qc_metrics,
                self.homer_make_ucsc_file,  #11
                self.macs2_callpeak,
                self.homer_annotate_peaks,
                self.homer_find_motifs_genome,
                self.annotation_graphs,
                self.run_spp,
                self.differential_binding, #17
                self.ihec_metrics,
                self.multiqc_report,
                self.cram_output,
                self.gatk_haplotype_caller,
                self.merge_and_call_individual_gvcf #22
            ],
            [
                self.picard_sam_to_fastq,
                self.trimmomatic,
                self.merge_trimmomatic_stats,
                self.mapping_bwa_mem_sambamba,
                self.sambamba_merge_bam_files,
                self.sambamba_mark_duplicates,
                self.sambamba_view_filter,
                self.metrics,
                self.homer_make_tag_directory,
                self.qc_metrics,
                self.homer_make_ucsc_file,
                self.macs2_atacseq_callpeak,
                self.homer_annotate_peaks,
                self.homer_find_motifs_genome,
                self.annotation_graphs,
                self.run_spp,
                self.differential_binding,
                self.ihec_metrics,
                self.multiqc_report,
                self.cram_output,
                self.gatk_haplotype_caller,
                self.merge_and_call_individual_gvcf
            ]
        ]


if __name__ == '__main__':

    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        ChipSeq(protocol=['chipseq', 'atacseq'])
