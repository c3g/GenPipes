#!/usr/bin/env python

# Python Standard Modules
import argparse
import logging
import os
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bio.design import *
from bio.readset import *

from pipelines.illumina import illumina

log = logging.getLogger(__name__)

# Those functions could be moved in a separate 'trinity' module in the future if they are reused elsewhere
def insilico_read_normalization(
    left_or_single_reads,
    right_reads,
    sequence_type,
    jellyfish_memory,
    output_directory=None,
    cpu=None
    ):

    normalization_stats_file = "normalization.stats.tsv"
    output_files = ["left.norm." + sequence_type, "right.norm." + sequence_type] if right_reads else ["single.norm." + sequence_type]
    if output_directory:
        output_files = [os.path.join(output_directory, output_file) for output_file in output_files]
        normalization_stats_file = os.path.join(output_directory, normalization_stats_file)

    output_files.append(normalization_stats_file)

    job = Job(
        left_or_single_reads + right_reads,
        output_files,
        [['insilico_read_normalization', 'module_perl'], ['insilico_read_normalization', 'module_trinity']]
    )

    job.command = """\
insilico_read_normalization.pl {other_options} \\
  --seqType {sequence_type} \\
  --JM {jellyfish_memory} \\
  --max_cov {maximum_coverage} \\
  {left_or_single_reads}{right_reads}{output_directory}{cpu}""".format(
        other_options=config.param('insilico_read_normalization', 'other_options', required=False),
        sequence_type=sequence_type,
        jellyfish_memory=jellyfish_memory,
        maximum_coverage=config.param('insilico_read_normalization', 'maximum_coverage', type="int"),
        left_or_single_reads=" \\\n  ".join(["--left " + read for read in left_or_single_reads]) if right_reads else " \\\n  ".join(["--single " + read for read in left_or_single_reads]),
        right_reads="".join([" \\\n  --right " + read for read in right_reads]) if right_reads else "",
        output_directory=" \\\n  --output " + output_directory if output_directory else "",
        cpu=" \\\n  --CPU " + str(cpu) if cpu else ""
    )

    if output_directory:
        job = concat_jobs([Job(command="mkdir -p " + output_directory), job])

    # Count normalized reads for stats
    job = concat_jobs([job, Job(command="""\
wc -l {output_file} | awk '{{print \\"# normalized {read_type} reads\t\\"\\$1 / 4}}' > {normalization_stats_file}""".format(
        output_file=output_files[0],
        read_type="paired" if right_reads else "single",
        normalization_stats_file=normalization_stats_file
    ))])

    return job


class RnaSeqDeNovoAssembly(illumina.Illumina):

    @property
    def contrasts(self):
        if not hasattr(self, "_contrasts"):
            self._contrasts = parse_design_file(self.args.design.name, self.samples)
        return self._contrasts

    def insilico_read_normalization_readsets(self):
        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            normalization_directory = os.path.join("insilico_read_normalization", readset.name)

            if readset.run_type == "PAIRED_END":
                left_or_single_reads = [trim_file_prefix + "pair1.fastq.gz"]
                right_reads = [trim_file_prefix + "pair2.fastq.gz"]
            elif readset.run_type == "SINGLE_END":
                left_or_single_reads = [trim_file_prefix + "single.fastq.gz"]
                right_reads = []
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            job = insilico_read_normalization(
                left_or_single_reads,
                right_reads,
                "fq",
                config.param('insilico_read_normalization_readsets', 'jellyfish_memory'),
                normalization_directory,
                config.param('insilico_read_normalization_readsets', 'cpu', required=False, type='int')
            )

            job.name = "insilico_read_normalization_readsets." + readset.name
            jobs.append(job)

        return jobs

    def insilico_read_normalization_all(self):
        normalization_directory = "insilico_read_normalization"
        left_or_single_reads = []
        right_reads = []

        for readset in self.readsets:
            if readset.run_type == "PAIRED_END":
                left_or_single_reads.append(os.path.join(normalization_directory, readset.name, "left.norm.fq"))
                right_reads.append(os.path.join(normalization_directory, readset.name, "right.norm.fq"))
            elif readset.run_type == "SINGLE_END":
                left_or_single_reads.append(trim_file_prefix + "single.norm.fq")
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

        job = insilico_read_normalization(
            left_or_single_reads,
            right_reads,
            "fq",
            config.param('insilico_read_normalization_all', 'jellyfish_memory'),
            os.path.join(normalization_directory, "all"),
            config.param('insilico_read_normalization_all', 'cpu', required=False, type='int')
        )

        job.name = "insilico_read_normalization_readsets." + readset.name
        return [job]

    def trinity(self):
        normalization_directory = os.path.join("insilico_read_normalization", "all")
        output_directory = "trinity_out_dir"
        trinity_fasta = os.path.join(output_directory, "Trinity.fasta")
        trinity_stats_prefix = os.path.join(output_directory, "Trinity.stats")

        if self.run_type == "PAIRED_END":
            left_reads = os.path.join(normalization_directory, "left.norm.fq")
            right_reads = os.path.join(normalization_directory, "right.norm.fq")
            input_files = [left_reads, right_reads]
            reads_option = "--left " + left_reads + " \\\n  --right " + right_reads
        elif self.run_type == "SINGLE_END":
            single_reads = os.path.join(normalization_directory, "single.norm.fq")
            input_files = [single_reads]
            reads_option = "--single " + single_reads

        trinity_job = Job(input_files, [trinity_fasta], [
            ['trinity', 'module_java'],
            ['trinity', 'module_trinity'],
            ['trinity', 'module_bowtie'],
            ['trinity', 'module_samtools']
        ])

        trinity_job.command = """\
Trinity {other_options} \\
  --JM {jellyfish_memory} \\
  --CPU {cpu} \\
  --bflyCPU {butterfly_cpu} \\
  {reads_option} \\
  --output {output_directory}""".format(
            other_options=config.param('trinity', 'other_options'),
            jellyfish_memory=config.param('trinity', 'jellyfish_memory'),
            cpu=config.param('trinity', 'cpu'),
            butterfly_cpu=config.param('trinity', 'butterfly_cpu'),
            reads_option=reads_option,
            output_directory=output_directory
        )

        return [concat_jobs([
            trinity_job,
            Job([trinity_fasta], [trinity_fasta + ".zip"], command="zip -j " + trinity_fasta + ".zip " + trinity_fasta),
            Job([trinity_fasta], [trinity_stats_prefix + ".csv", trinity_stats_prefix + ".jpg", trinity_stats_prefix + ".pdf"], [['trinity', 'module_R']], command="Rscript -e 'library(gqSeqUtils); dnaFastaStats(filename = \\\"" + trinity_fasta + "\\\", type = \\\"trinity\\\", output.prefix = \\\"" + trinity_stats_prefix + "\\\")'")
        ], name="trinity")]

    def exonerate_fastasplit(self):
        trinity_directory = "trinity_out_dir"
        trinity_fasta = os.path.join(trinity_directory, "Trinity.fasta")
        trinity_chunks_directory = os.path.join(trinity_directory, "Trinity.fasta_chunks")
        num_fasta_chunks = config.param('exonerate_fastasplit', 'num_fasta_chunks', type='posint')

        return [concat_jobs([
            Job(command="mkdir -p " + trinity_chunks_directory),
            Job(
                [trinity_fasta],
                [os.path.join(trinity_chunks_directory, "Trinity.fasta_chunk_{:07d}".format(i)) for i in range(num_fasta_chunks)],
                [['exonerate_fastasplit', 'module_exonerate']],
                command="fastasplit -f " + trinity_fasta + " -o " + trinity_chunks_directory + " -c " + str(num_fasta_chunks)
            )
        ], name="exonerate_fastasplit.Trinity.fasta")]

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.insilico_read_normalization_readsets,
            self.insilico_read_normalization_all,
            self.trinity,
            self.exonerate_fastasplit
        ]

    def __init__(self):
        # Add pipeline specific arguments
        self.argparser.add_argument("-d", "--design", help="design file", type=file)

        super(RnaSeqDeNovoAssembly, self).__init__()
        
RnaSeqDeNovoAssembly().submit_jobs()
