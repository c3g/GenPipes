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
import os
import sys
import logging
import collections
import re

# Append mugqic_pipelines directory to Python library path
import distutils.util
from typing import List, Any

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
import utils.utils
from core.config import global_config_parser, SanitycheckError, _raise
from core.job import Job, concat_jobs, pipe_jobs
from bfx.readset import parse_nanopore_readset_file
from pipelines import common

from bfx import artic
from bfx import covseq_tools
from bfx import guppy
from bfx import htslib
from bfx import kraken2
from bfx import minimap2
from bfx import ncovtools
from bfx import pycoqc
from bfx import quast
from bfx import sambamba
from bfx import samtools
from bfx import snpeff
from bfx import wub

from bfx import bash_cmd as bash

log = logging.getLogger(__name__)

class NanoporeCoVSeq(common.Nanopore):
    """
    Nanopore CoVSeq Pipeline
    ==============


    For information on the structure and contents of the Nanopore readset file, please consult [here](
    https://bitbucket.org/mugqic/genpipes/src/master/#markdown-header-nanopore).
    """

    def __init__(self, *args, protocol=None, **kwargs ):
        self._protocol = protocol
        super(NanoporeCoVSeq, self).__init__(*args, **kwargs)

    @classmethod
    def argparser(cls, argparser):
        super().argparser(argparser)
        cls._argparser.add_argument("-t", "--type",
                                    help="Type of CoVSeQ analysis,basecalling on/off (default without basecalling)",
                                    choices=["default", "basecalling"], default="default", dest='protocol')

        return cls._argparser

    @property
    def samples(self):
        if self._samples is None:
            self._samples = list(collections.OrderedDict.fromkeys([readset.sample for readset in self.readsets]))
            reads_fast5_dir = []
            for sample in self.samples:
                reads_fast5_dir.append(sample.fast5_files)
                if not sample.fast5_files:
                    _raise(SanitycheckError("Error: FAST5 file not available for readset \"" + sample.name + "\"!"))
                    exit()

            reads_fast5_dir = list(set(reads_fast5_dir))
            if len(reads_fast5_dir) > 1:
                _raise(SanitycheckError("Error: FAST5 files are organized in more than 1 folder: \"" + "\n\t- ".join(reads_fast5_dir) + "\"!"))
                exit()
        return self._samples

    @property
    def run_name(self):
        if not hasattr(self, "_run_name"):
            self._run_name = global_config_parser.param('DEFAULT', 'run_name', required=True)
        return self._run_name

    def guppy_basecall(self):
        """
        Use guppy to perform basecalling on raw FAST5 files

        """
        jobs = []

        reads_fast5_dir = []
        transfer = bool(distutils.util.strtobool(global_config_parser.param('guppy_basecall', 'transfer_to_tmp')))

        for sample in self.samples:
            reads_fast5_dir.append(sample.fast5_files)

        reads_fast5_dir = list(set(reads_fast5_dir))
        fastq_directory = os.path.join("basecall")

        jobs.append(
            concat_jobs([
                guppy.guppy_basecalling(
                    reads_fast5_dir,
                    fastq_directory,
                    transfer
                ),
                Job(
                    input_files=[],
                    output_files=[os.path.join(fastq_directory, "logs.tar.gz")],
                    command="""tar --remove-files -acvf {fastq_directory}/logs.tar.gz {fastq_directory}/*.log""".format(
                        fastq_directory=fastq_directory
                    )
                )
            ],
                name="guppy_basecall"
            )
        )

        fastq_directory = os.path.join("basecall")

        job = guppy.guppy_basecalling(reads_fast5_dir, fastq_directory)
        job.name = "guppy_basecall"

        return jobs


    def guppy_demultiplex(self):
        """
        Use guppy to perform demultiplexing on raw FASTQ read files

        """
        jobs = []
        fastq_directory = os.path.join("basecall", "pass")
        sequencing_summary = os.path.join("basecall", "sequencing_summary.txt")

        demux_fastq_directory = os.path.join("demultiplex")

        transfer = bool(distutils.util.strtobool(global_config_parser.param('guppy_demultiplex', 'transfer_to_tmp')))

        demux_barcode_dir = []
        for sample in self.samples:
            if not sample.barcode:
                _raise(SanitycheckError("Error: Sample \"" + sample.name + "\" doesn't have a barcode associated!"))
            demux_barcode_dir.append(os.path.join(demux_fastq_directory, sample.barcode))

        jobs.append(
            concat_jobs([
                guppy.guppy_demultiplex(
                    fastq_directory,
                    sequencing_summary,
                    demux_fastq_directory,
                    demux_barcode_dir,
                    transfer
                ),
                Job(
                    input_files=[],
                    output_files=[os.path.join(demux_fastq_directory, "logs.tar.gz")],
                    command="""tar --remove-files -acvf {demux_fastq_directory}/logs.tar.gz {demux_fastq_directory}/*.log""".format(
                        demux_fastq_directory=demux_fastq_directory
                    )
                )
            ],
                name="guppy_demultiplex"
            )
        )

        return jobs

    def pycoqc(self):
        """
        Use pycoQC to produce an interactive quality report based on the summary file and
        alignment outputs.
        """
        jobs = []

        fastq_directory = os.path.join("basecall")
        sequencing_summary = os.path.join(fastq_directory, "sequencing_summary.txt")
        pycoqc_directory = os.path.join("report", "pycoQC")
        # run_name = self.run_name

        run_names = []
        for sample in self.samples:
            run_names.append(sample.run)

        for run_name in list(set(run_names)):
            jobs.append(
                concat_jobs([
                    bash.mkdir(pycoqc_directory),
                    pycoqc.pycoqc(
                        run_name,
                        sequencing_summary,
                        pycoqc_directory
                    )
                ],
                    name="pycoqc"
                )
            )

        # job.name = "pycoqc_report"

        return jobs

    def host_reads_removal_dependency(self):
        """
        Runs minimap2 on a hybrid genome to remove potential host reads
        """

        jobs = []

        demux_fastq_directory = os.path.join("demultiplex")

        for sample in self.samples:
            host_removal_directory = os.path.join("host_removal", sample.name)
            reads_fastq_dir = os.path.join(demux_fastq_directory, sample.barcode)
            sample_bam = os.path.join(host_removal_directory, sample.name + ".hybrid.sorted.bam")
            sample_bam_host_removed_sorted = os.path.join(host_removal_directory, sample.name + ".host_removed.sorted.bam")
            sample_bam_host_removed_sorted_index = os.path.join(host_removal_directory, sample.name + ".host_removed.sorted.bam.bai")
            output_fq = os.path.join(host_removal_directory, sample.name + ".host_removed.fastq.gz")
            jobs.append(
                concat_jobs([
                    bash.mkdir(host_removal_directory),
                    pipe_jobs([
                        minimap2.minimap2_ont(
                            reads_fastq_dir,
                            read_group="'@RG" + \
                                       "\\tID:" + sample.name + \
                                       "\\tSM:" + sample.name + \
                                       "\\tLB:" + (sample.library if sample.library else sample.name) + \
                                       ("\\tPU:run" + sample.run if sample.run else "") + \
                                       "\\tPL:Nanopore" + \
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
                            sample_bam,
                            tmp_dir=global_config_parser.param('host_reads_removal', 'tmp_dir', required=True),
                            other_options=global_config_parser.param('host_reads_removal', 'sambamba_sort_other_options',
                                                                     required=False)
                        )
                    ]),
                    sambamba.view(
                        sample_bam,
                        sample_bam_host_removed_sorted,
                        options=global_config_parser.param('host_reads_removal', 'sambamba_view_other_options')
                    ),
                    sambamba.index(
                        sample_bam_host_removed_sorted,
                        sample_bam_host_removed_sorted_index,
                        other_options=global_config_parser.param('host_reads_removal', 'sambamba_index_other_options', required=False)
                    ),
                    samtools.bam2fq(
                        input_bam=sample_bam_host_removed_sorted,
                        output_pair1=None,
                        output_pair2=None,
                        output_other=output_fq,
                        output_single=None,
                        ini_section='host_reads_removal'
                    )
                ],
                    name="host_reads_removal." + sample.name,
                    samples=[sample]
                )
            )

        return jobs

    def host_reads_removal(self):
        """
        Runs minimap2 on a hybrid genome to remove potential host reads
        """

        jobs = []

        # demux_fastq_directory = os.path.join("demultiplex")

        for sample in self.samples:
            if sample.fastq_files:
                demux_fastq_directory = sample.fastq_files
            else:
                _raise(SanitycheckError("Error: FASTQ files not available for sample \"" + sample.name + "\"!"))

            host_removal_directory = os.path.join("host_removal", sample.name)
            reads_fastq_dir = os.path.join(demux_fastq_directory, sample.barcode)
            sample_bam = os.path.join(host_removal_directory, sample.name + ".hybrid.sorted.bam")
            sample_bam_host_removed_sorted = os.path.join(host_removal_directory, sample.name + ".host_removed.sorted.bam")
            sample_bam_host_removed_sorted_index = os.path.join(host_removal_directory, sample.name + ".host_removed.sorted.bam.bai")
            output_fq = os.path.join(host_removal_directory, sample.name + ".host_removed.fastq.gz")
            jobs.append(
                concat_jobs([
                    bash.mkdir(host_removal_directory),
                    pipe_jobs([
                        minimap2.minimap2_ont(
                            reads_fastq_dir,
                            read_group="'@RG" + \
                                       "\\tID:" + sample.name + \
                                       "\\tSM:" + sample.name + \
                                       "\\tLB:" + (sample.library if sample.library else sample.name) + \
                                       ("\\tPU:run" + sample.run if sample.run else "") + \
                                       "\\tPL:Nanopore" + \
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
                            sample_bam,
                            tmp_dir=global_config_parser.param('host_reads_removal', 'tmp_dir', required=True),
                            other_options=global_config_parser.param('host_reads_removal', 'sambamba_sort_other_options',
                                                                     required=False)
                        )
                    ]),
                    sambamba.view(
                        sample_bam,
                        sample_bam_host_removed_sorted,
                        options=global_config_parser.param('host_reads_removal', 'sambamba_view_other_options')
                    ),
                    sambamba.index(
                        sample_bam_host_removed_sorted,
                        sample_bam_host_removed_sorted_index,
                        other_options=global_config_parser.param('host_reads_removal', 'sambamba_index_other_options', required=False)
                    ),
                    samtools.bam2fq(
                        input_bam=sample_bam_host_removed_sorted,
                        output_pair1=None,
                        output_pair2=None,
                        output_other=output_fq,
                        output_single=None,
                        ini_section='host_reads_removal'
                    )
                ],
                    name="host_reads_removal." + sample.name,
                    samples=[sample]
                )
            )

        return jobs

    def kraken_analysis(self):
        """
        kraken
        """

        jobs = []
        for sample in self.samples:
            host_removal_directory = os.path.join("host_removal", sample.name)
            kraken_directory = os.path.join("metrics", "dna", sample.name, "kraken_metrics")
            kraken_out_prefix = os.path.join(kraken_directory, sample.name)

            [fastq1] = [os.path.join(host_removal_directory, sample.name + ".host_removed.fastq.gz")]
            fastq2 = None
            unclassified_output = [kraken_out_prefix + ".unclassified_sequences.fastq"]
            classified_output = [kraken_out_prefix + ".classified_sequences.fastq"]

            jobs.append(
                concat_jobs([
                    bash.mkdir(kraken_directory),
                    kraken2.kraken2(
                        fastq1,
                        fastq2,
                        kraken_out_prefix,
                        other_options=global_config_parser.param('kraken_analysis', 'kraken2_other_options'),
                        nthread=global_config_parser.param('kraken_analysis', 'kraken2_threads'),
                        database=global_config_parser.param('kraken_analysis', 'kraken2_database')
                    ),
                    Job(
                        input_files=unclassified_output + classified_output,
                        output_files=[s + ".gz" for s in unclassified_output + classified_output],
                        module_entries=[
                            ['pigz', 'module_pigz']
                        ],
                        command="""pigz -k -f -p {nthreads} {input_files}""".format(
                            input_files=" ".join(unclassified_output + classified_output),
                            nthreads=global_config_parser.param('kraken_analysis', 'pigz_threads')
                        )
                    )
                ],
                    name="kraken_analysis." + sample.name,
                    removable_files=unclassified_output + classified_output
                )
            )

        return jobs


    def artic_nanopolish(self):
        """
        Runs artic nanopolish pipeline on all samples.
        """

        jobs = []

        for sample in self.samples:

            artic_nanopolish_directory = os.path.join("artic_nanopolish", sample.name)
            reads_fastq_dir = os.path.join("host_removal", sample.name)
            variant_directory = os.path.join("variant", sample.name)
            variant_filename = sample.name + ".pass.vcf.gz"
            variant = os.path.join(variant_directory, variant_filename)
            variant_artic = os.path.join(artic_nanopolish_directory, variant_filename)
            variant_link = os.path.join("..", "..", artic_nanopolish_directory, variant_filename)
            variant_index_filename = sample.name + ".pass.vcf.gz.tbi"
            variant_index = os.path.join(variant_directory, variant_index_filename)
            variant_index_artic = os.path.join(artic_nanopolish_directory, variant_index_filename)
            variant_index_link = os.path.join("..", "..", artic_nanopolish_directory, variant_index_filename)
            consensus_directory = os.path.join("consensus", sample.name)
            consensus_filename = sample.name + ".consensus.fasta"
            consensus = os.path.join(consensus_directory, consensus_filename)
            consensus_artic = os.path.join(artic_nanopolish_directory, consensus_filename)
            consensus_link = os.path.join("..", "..", artic_nanopolish_directory, consensus_filename)

            alignment_directory = os.path.join("alignment", sample.name)
            # bam before primer trimming
            raw_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            raw_bam_artic = os.path.join(artic_nanopolish_directory, sample.name + ".trimmed.rg.sorted.bam")
            raw_bam_link = os.path.join("..", "..", artic_nanopolish_directory, sample.name + ".trimmed.rg.sorted.bam")
            raw_bam_index = os.path.join(alignment_directory, sample.name + ".sorted.bam.bai")
            raw_bam_index_artic = os.path.join(artic_nanopolish_directory, sample.name + ".trimmed.rg.sorted.bam.bai")
            raw_bam_index_link = os.path.join("..", "..", artic_nanopolish_directory,
                                              sample.name + ".trimmed.rg.sorted.bam.bai")
            # bam after primer trimming
            primer_trimmed_bam = os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.bam")
            primer_trimmed_bam_artic = os.path.join(artic_nanopolish_directory,
                                                    sample.name + ".primertrimmed.rg.sorted.bam")
            primer_trimmed_bam_index = os.path.join(alignment_directory,
                                                    sample.name + ".sorted.filtered.primerTrim.bam.bai")

            if sample.summary_file:
                sequencing_summary = sample.summary_file
            else:
                sequencing_summary = os.path.join("basecall", "sequencing_summary.txt")

            jobs.append(
                concat_jobs([
                    bash.mkdir(artic_nanopolish_directory),
                    bash.chgdir(artic_nanopolish_directory),
                    artic.nanopolish_ont(
                        reads_fastq_dir,
                        sample.run,
                        sample.name,
                        sample.fast5_files,
                        sequencing_summary,
                        artic_nanopolish_directory,
                        ini_section="artic_nanopolish"
                    ),
                    bash.chgdir(self.output_dir),
                    bash.mkdir(consensus_directory),
                    Job(
                        input_files=[consensus_artic],
                        output_files=[consensus],
                        command="""ln -sf {consensus_link} {consensus}""".format(
                            consensus_link=consensus_link,
                            consensus=consensus
                        )
                    ),
                    bash.mkdir(variant_directory),
                    Job(
                        input_files=[variant_artic],
                        output_files=[variant],
                        command="""ln -sf {variant_link} {variant}""".format(
                            variant_link=variant_link,
                            variant=variant
                        )
                    ),
                    Job(
                        input_files=[variant_index_artic],
                        output_files=[variant_index],
                        command="""ln -sf {variant_index_link} {variant_index}""".format(
                            variant_index_link=variant_index_link,
                            variant_index=variant_index
                        )
                    ),
                    bash.mkdir(alignment_directory),
                    Job(
                        input_files=[raw_bam_artic],
                        output_files=[raw_bam],
                        command="""ln -sf {raw_bam_link} {raw_bam}""".format(
                            raw_bam_link=raw_bam_link,
                            raw_bam=raw_bam
                        )
                    ),
                    Job(
                        input_files=[raw_bam_index_artic],
                        output_files=[raw_bam_index],
                        command="""ln -sf {raw_bam_index_link} {raw_bam_index}""".format(
                            raw_bam_index_link=raw_bam_index_link,
                            raw_bam_index=raw_bam_index
                        )
                    ),
                    pipe_jobs([
                        sambamba.view(
                            primer_trimmed_bam_artic,
                            None,
                            "-f bam -F \"not supplementary and not secondary_alignment\""
                        ),
                        sambamba.sort(
                            "/dev/stdin",
                            primer_trimmed_bam,
                            tmp_dir=global_config_parser.param('artic_nanopolish', 'tmp_dir', required=True)
                        )
                    ]),
                    sambamba.index(
                        primer_trimmed_bam,
                        primer_trimmed_bam_index
                    )
                ],
                    name="artic_nanopolish." + sample.name
                ),
            )

        return jobs

    def wub_metrics(self):
        """
        Generate WUB metrics on bam file
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            metrics_dir = os.path.join("metrics", "dna", sample.name, "wub_metrics")
            primer_trimmed_bam = os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.bam")
            pickle_out = os.path.join(metrics_dir, sample.name + ".bam_alignment_qc.pk")

            jobs.append(
                concat_jobs([
                    bash.mkdir(metrics_dir),
                    wub.bam_alignment_qc(
                        primer_trimmed_bam,
                        pickle_out,
                        ini_section="wub_metrics"
                    ),
                ],
                    name="wub_metrics." + sample.name,
                    samples=[sample]
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

            input_vcf = os.path.join(variant_directory, sample.name + ".pass.vcf.gz")
            output_vcf = os.path.join(variant_directory,
                                      re.sub("\.vcf.gz$", ".annotate.vcf", os.path.basename(input_vcf)))

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
                    ),
                ],
                    name="snpeff_annotate." + sample.name,
                    samples=[sample]
                )
            )

        return jobs

    def covseq_metrics(self):
        """
        """

        jobs = []

        for sample in self.samples:
            consensus_directory = os.path.join("consensus", sample.name)
            consensus = os.path.join(consensus_directory, sample.name + ".consensus.fasta")
            metrics_directory = os.path.join("metrics", "dna", sample.name, "general_metrics")
            output = os.path.join(metrics_directory, sample.name + ".metrics.csv")
            fq_stats = os.path.join(metrics_directory, sample.name + ".fastq.stats")
            artic_nanopolish_directory = os.path.join("artic_nanopolish", sample.name)
            pass_fq = os.path.join(artic_nanopolish_directory, sample.run + "_" + sample.name + ".fastq")
            wub_metrics_dir = os.path.join("metrics", "dna", sample.name, "wub_metrics")
            pickle = os.path.join(wub_metrics_dir, sample.name + ".bam_alignment_qc.pk")
            if sample.summary_file:
                sequencing_summary = sample.summary_file
            else:
                sequencing_summary = os.path.join("basecall", "sequencing_summary.txt")
            jobs.append(
                concat_jobs([
                    bash.mkdir(metrics_directory),
                    Job(
                        input_files=[sequencing_summary, pass_fq],
                        output_files=[fq_stats],
                        command="""\\
echo "raw_reads" $(grep -c {barcode} {sequencing_summary}) > {fq_stats} && \\
echo "pass_reads" $(grep -c "^@" {pass_fq}) >> {fq_stats} """.format(
                            barcode=sample.barcode,
                            sequencing_summary=sequencing_summary,
                            sample=sample.name,
                            pass_fq=pass_fq,
                            fq_stats=fq_stats
                        )
                    ),
                    covseq_tools.covid_collect_nanopore_metrics(
                        sample.name,
                        consensus,
                        fq_stats,
                        pickle,
                        output,
                        ini_section='covseq_metrics'
                    )
                ],
                    name="covseq_metrics." + sample.name,
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
            consensus = os.path.join(consensus_directory, sample.name + ".consensus.fasta")

            jobs.append(
                concat_jobs([
                    bash.mkdir(consensus_directory),
                    quast.quast(
                        consensus,
                        output_dir
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
        for sample in self.samples:
            consensus_directory = os.path.join("consensus", sample.name)
            input_fa = os.path.join(consensus_directory, sample.name + ".consensus.fasta")
            quast_directory = os.path.join("metrics", "dna", sample.name, "quast_metrics")
            quast_html = os.path.join(quast_directory, "report.html")
            quast_tsv = os.path.join(quast_directory, "report.tsv")

            output_fa = os.path.join(consensus_directory, sample.name + ".consensus.fasta")
            output_status_fa = os.path.join(consensus_directory,
                                            """{sample_name}.consensus.{technology}.{status}.fasta""".format(
                                                sample_name=sample.name,
                                                technology=global_config_parser.param('rename_consensus_header',
                                                                        'sequencing_technology', required=False),
                                                status="${STATUS}"))

            variant_directory = os.path.join("variant", sample.name)
            input_vcf = os.path.join(variant_directory, sample.name + ".pass.vcf.gz")
            annotated_vcf = os.path.join(variant_directory,
                                         re.sub("\.vcf.gz$", ".annotate.vcf", os.path.basename(input_vcf)))

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(output_fa)),
                    Job(
                        input_files=[quast_tsv, quast_html],
                        output_files=[],
                        command="""\\
cons_len=`grep -oP "Total length \(>= 0 bp\)\\t\K.*?(?=$)" {quast_tsv}`
N_count=`grep -oP "# N's\\",\\"quality\\":\\"Less is better\\",\\"values\\":\[\K.*?(?=])" {quast_html}`
cons_perc_N=`echo "scale=2; 100*$N_count/$cons_len" | bc -l`
frameshift=`if grep -q "frameshift_variant" {annotated_vcf}; then echo "FLAG"; fi`
STATUS=`awk -v frameshift=$frameshift -v cons_perc_N=$cons_perc_N 'BEGIN {{ if (cons_perc_N < 5 && frameshift != "FLAG") {{print "pass"}} else if ((cons_perc_N >= 5 && cons_perc_N <= 10 ) || frameshift == "FLAG") {{print "flag"}} else if (cons_perc_N > 10) {{print "rej"}} }}'`
export STATUS""".format(
                            quast_html=quast_html,
                            quast_tsv=quast_tsv,
                            annotated_vcf=annotated_vcf
                        )
                    ),
                    Job(
                        input_files=[input_fa],
                        output_files=[output_status_fa],
                        command="""\\
awk '/^>/{{print ">{country}/{province}-{sample}/{year} seq_method:{seq_method}|assemb_method:{assemb_method}|snv_call_method:{snv_call_method}"; next}}{{print}}' < {input_fa} > {output_status_fa}""".format(
                            country=global_config_parser.param('rename_consensus_header', 'country', required=False),
                            province=global_config_parser.param('rename_consensus_header', 'province', required=False),
                            year=global_config_parser.param('rename_consensus_header', 'year', required=False),
                            seq_method=global_config_parser.param('rename_consensus_header', 'seq_method', required=False),
                            assemb_method=global_config_parser.param('rename_consensus_header', 'assemb_method', required=False),
                            snv_call_method=global_config_parser.param('rename_consensus_header', 'snv_call_method', required=False),
                            sample=sample.name,
                            input_fa=input_fa,
                            output_status_fa=output_status_fa
                        )
                    )
                ],
                    name="rename_consensus_header." + sample.name,
                    samples=[sample]
                )
            )

        return jobs

    def prepare_report(self):

        jobs = []

        readset_file = os.path.relpath(self.readsets_file.name, self.output_dir)
        readset_file_report = "report.readset.tsv"

        software_version = os.path.join("report", "software_versions.csv")
        run_metadata = os.path.join("report", "run_metadata.csv")

        ncovtools_directory = os.path.join("report", "ncov_tools")
        metadata = os.path.join(ncovtools_directory, "metadata.tsv")
        ncovtools_data_directory = os.path.join(ncovtools_directory, "data")
        ncovtools_config = os.path.join(ncovtools_directory, "config.yaml")

        modules = []
        # Retrieve all unique module version values in config files
        # assuming that all module key names start with "module_"
        for section in global_config_parser.sections():
            for name, value in global_config_parser.items(section):
                if re.search("^module_", name) and value not in modules:
                    modules.append(value)

        job = concat_jobs([
            bash.mkdir(ncovtools_data_directory),
            Job(
                input_files=[],
                output_files=[readset_file_report, metadata],
                command="""\\
head -n 1 {readset_file} > {readset_file_report} && \\
echo -e "sample\\tct\\tdate" > {metadata}""".format(
                    readset_file=readset_file,
                    readset_file_report=readset_file_report,
                    metadata=metadata
                )
            ),
        ])

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            raw_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            primer_trimmed_bam = os.path.join(alignment_directory, sample.name + ".sorted.filtered.primerTrim.bam")
            consensus_directory = os.path.join("consensus", sample.name)
            consensus = os.path.join(consensus_directory, sample.name + ".consensus.fasta")
            variant_directory = os.path.join("variant", sample.name)
            variants = os.path.join(variant_directory, sample.name + ".pass.vcf.gz")
            variants_index = os.path.join(variant_directory, sample.name + ".pass.vcf.gz.tbi")

            output_raw_bam = os.path.join(ncovtools_data_directory, os.path.basename(raw_bam))
            output_primer_trimmed_bam = os.path.join(ncovtools_data_directory, os.path.basename(primer_trimmed_bam))
            output_consensus = os.path.join(ncovtools_data_directory, os.path.basename(consensus))
            output_variants = os.path.join(ncovtools_data_directory, os.path.basename(variants))
            output_variants_index = os.path.join(ncovtools_data_directory, os.path.basename(variants_index))

            job = concat_jobs([
                job,
                Job(
                    input_files=[raw_bam, primer_trimmed_bam, consensus, variants],
                    output_files=[output_raw_bam, output_primer_trimmed_bam, output_consensus, output_variants],
                    command="""\\
echo "Linking files for ncov_tools for sample {sample_name}..." && \\
if [ "$(ls -1 {raw_bam})" != "" ] && [ "$(ls -1 {primer_trimmed_bam})" != "" ] && [ "$(ls -1 {consensus})" != "" ] && [ "$(ls -1 {variants})" != "" ];
  then
    ln -fs $(pwd -P )/$(ls -1 {raw_bam}) {output_raw_bam} && \\
    ln -fs $(pwd -P )/$(ls -1 {primer_trimmed_bam}) {output_primer_trimmed_bam} && \\
    ln -fs $(pwd -P )/$(ls -1 {consensus}) {output_consensus} && \\
    ln -fs $(pwd -P )/$(ls -1 {variants}) {output_variants} && \\
    ln -fs $(pwd -P )/$(ls -1 {variants_index}) {output_variants_index} && \\
    grep {sample_name} {readset_file} >> {readset_file_report} && \\
    echo -e "{sample_name}\\tNA\\tNA" >> {metadata}
fi""".format(
                        readset_file=readset_file,
                        readset_file_report=readset_file_report,
                        raw_bam=raw_bam,
                        primer_trimmed_bam=primer_trimmed_bam,
                        consensus=consensus,
                        variants=variants,
                        variants_index=variants_index,
                        output_raw_bam=output_raw_bam,
                        output_primer_trimmed_bam=output_primer_trimmed_bam,
                        output_consensus=output_consensus,
                        output_variants=output_variants,
                        output_variants_index=output_variants_index,
                        sample_name=sample.name,
                        metadata=metadata
                    )
                )
            ],
                samples=[sample]
            )
        jobs.append(
            concat_jobs([
                job,
                ncovtools.run_ncovtools(
                    output_raw_bam,
                    output_primer_trimmed_bam,
                    output_consensus,
                    output_variants,
                    readset_file,
                    metadata,
                    ncovtools_directory,
                    ncovtools_config,
                    self.output_dir,
                )
            ],
                name="prepare_report." + self.run_name
            )
        )

        return jobs

    @property
    def step_list(self):
        return self.protocols()[self._protocol]

    def protocols(self):
        return {
            'default':
            [
                self.host_reads_removal,
                self.kraken_analysis,
                self.artic_nanopolish,
                self.wub_metrics,
                self.covseq_metrics,
                self.snpeff_annotate,
                self.quast_consensus_metrics,
                self.rename_consensus_header,
                self.prepare_report],
            'basecalling': [
                self.guppy_basecall,
                self.guppy_demultiplex,
                self.pycoqc,
                self.host_reads_removal_dependency,
                self.kraken_analysis,
                self.artic_nanopolish,
                self.wub_metrics,
                self.covseq_metrics,
                self.snpeff_annotate,
                self.quast_consensus_metrics,
                self.rename_consensus_header,
                self.prepare_report],
        }


def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    # Check if Genpipes must be ran inside a container
    utils.container_wrapper_argparse(__file__, argv)
    # Build help
    epilog = NanoporeCoVSeq.process_help(argv)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve', epilog=epilog)

    # populate the parser
    parser = NanoporeCoVSeq.argparser(parser)

    parsed_args = parser.parse_args(argv)

    sanity_check = parsed_args.sanity_check
    loglevel = parsed_args.log
    utils.set_logger(loglevel, sanity_check=sanity_check)

    # Pipeline config
    config_files = parsed_args.config

    # Common Pipeline options
    genpipes_file = parsed_args.genpipes_file
    container = parsed_args.container
    clean = parsed_args.clean
    report = parsed_args.report
    no_json = parsed_args.no_json
    force = parsed_args.force
    job_scheduler = parsed_args.job_scheduler
    output_dir = parsed_args.output_dir
    steps = parsed_args.steps
    readset_file = parsed_args.readsets_file
    design_file = parsed_args.design_file
    protocol = parsed_args.protocol

    pipeline = NanoporeCoVSeq(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file,
                              clean=clean, report=report, force=force, job_scheduler=job_scheduler, output_dir=output_dir,
                              design_file=design_file, no_json=no_json, container=container,
                              protocol=protocol)

    pipeline.submit_jobs()

if __name__ == '__main__':
    main()
