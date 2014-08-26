#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def cuffcompare(gtf_files, output_prefix, gtf_list):
    job = Job(
        gtf_files,
        [output_prefix + ".combined.gtf", output_prefix + ".TranscriptList.tsv"],
        [['cuffcompare', 'module_cufflinks'], ['cuffcompare', 'module_tools']]
    )

    job.command = """\
mkdir -p {output_directory} && \\
cuffcompare -T \\
  -o {output_prefix} \\
  -r {reference_gff} \\
  -R {reference_fasta} \\
  {gtf_files} && \\
formatDenovoCombinedGTF.py \\
  -c {output_prefix}.combined.gtf \\
  -t {output_prefix}.tracking \\
  -s {gtf_list} \\
  -o {output_prefix}.TranscriptList.tsv""".format(
        output_directory=os.path.dirname(output_prefix),
        output_prefix=output_prefix,
        reference_gff=config.param('cuffcompare', 'gtf', type='filepath'),
        reference_fasta=config.param('cuffcompare', 'genome_fasta', type='filepath'),
        gtf_files=" \\\n  ".join(gtf_files),
        gtf_list=gtf_list
    )

    return job

def cuffdiff(sample_replicate_group_bams, gtf, output_directory):

    # sample_replicate_group_bams is a list of lists of replicates per sample
    # Flatten this list to set job input files
    input_bams = []
    for sample_replicate_bams in sample_replicate_group_bams:
        input_bams.extend(sample_replicate_bams)

    job = Job(
        input_bams + [gtf],
        [os.path.join(output_directory, "isoform_exp.diff")],
        [['cuffcompare', 'module_cufflinks']]
    )

    job.command = """\
mkdir -p {output_directory} && \\
cuffdiff {other_options} \\
  --frag-bias-correct {genome_fasta} \\
  --library-type {library_type} \\
  --output-dir {output_directory} \\
  --num-threads {num_threads} \\
  {gtf} \\
  {input_bams}""".format(
        other_options=config.param('cuffdiff', 'other_options'),
        genome_fasta=config.param('cuffdiff', 'genome_fasta', type='filepath'),
        library_type=config.param('cuffdiff', 'library_type'),
        output_directory=output_directory,
        num_threads=config.param('cuffdiff', 'threads', type='posint'),
        gtf=gtf,
        # Join replicate bams per sample with a "," then join all sample replicate groups with a " "
        input_bams=" \\\n  ".join([",".join(sample_replicate_bams) for sample_replicate_bams in sample_replicate_group_bams])
    )

    return job

def cufflinks(input_bam, output_directory, gtf=None):

    job = Job(
        [input_bam],
        [os.path.join(output_directory, "transcripts.gtf")],
        [['cufflinks', 'module_cufflinks']]
    )

    job.command = """\
mkdir -p {output_directory} && \\
cufflinks -q {other_options}{gtf} \\
  --max-bundle-frags {max_bundle_frags} \\
  --library-type {library_type} \\
  --output-dir {output_directory} \\
  --num-threads {num_threads} \\
  {input_bam}""".format(
        other_options=config.param('cufflinks', 'other_options', required=False),
        gtf=" \\\n  --GTF " + gtf if gtf else "",
        max_bundle_frags=config.param('cufflinks', 'max_bundle_frags', type='int'),
        library_type=config.param('cufflinks', 'library_type'),
        output_directory=output_directory,
        num_threads=config.param('cufflinks', 'threads', type='posint'),
        input_bam=input_bam
    )

    return job
