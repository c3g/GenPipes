#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def dna_sample_metrics(input_directory, output, experiment_type="unknown"):
    job = Job([input_directory], [output], [['dna_sample_metrics', 'moduleVersion.R'], ['dna_sample_metrics', 'moduleVersion.tools']])

    job.command = \
"""Rscript \$R_TOOLS/DNAsampleMetrics.R \\
  {input_directory} \\
  {output} \\
  {experiment_type}""".format(
        input_directory=input_directory,
        output=output,
        experiment_type=experiment_type
    )

    return job

def merge_trimmomatic_stats(input_pattern, input_directory, output, type):
    job = Job([], [output], [['merge_trimmomatic_stats', 'moduleVersion.R'], ['merge_trimmomatic_stats', 'moduleVersion.tools']])

    job.command = \
"""mkdir -p {output_directory} && \\
Rscript \$R_TOOLS/mergeTrimmomaticStat.R \\
  {input_pattern} \\
  {input_directory} \\
  {output} \\
  {type}""".format(
        output_directory=os.path.dirname(output),
        input_pattern=input_pattern,
        input_directory=input_directory,
        output=output,
        type=type
    )

    return job

def rnaseqc(sample_file, output_directory, run_type=None):
    job = Job([sample_file], [os.path.join(output_directory, "index.html")], [['rnaseqc', 'moduleVersion.java'], ['rnaseqc', 'moduleVersion.bwa'], ['rnaseqc', 'moduleVersion.rnaseqc']])

    job.command = \
"""java -Djava.io.tmpdir={tmp_dir} {extra_java_flags} -Xmx{ram} -jar \$RNASEQC_JAR \\
  -BWArRNA {reference_ribosomal_rna_fasta} \\
  -n {number_top_transcripts} \\
  -o {output_directory} \\
  -r {reference_genome_fasta} \\
  -s {sample_file} \\
  -t {gtf_file}{single_end}""".format(
        tmp_dir=config.param('rnaseqc', 'tmpDir'),
        extra_java_flags=config.param('rnaseqc', 'extraJavaFlags'),
        ram=config.param('rnaseqc', 'ram'),
        reference_ribosomal_rna_fasta=config.param('rnaseqc', 'ribosomalFasta', type='filepath'),
        number_top_transcripts=config.param('rnaseqc', 'topTranscript', type='int'),
        output_directory=output_directory,
        reference_genome_fasta=config.param('rnaseqc', 'referenceFasta', type='filepath'),
        sample_file=sample_file,
        gtf_file=config.param('rnaseqc', 'referenceGtf', type='filepath'),
        single_end=" \\\n  -singleEnd" if run_type == "SINGLE_END" else ""
    )

    return job

def snv_graph_metrics(list, output_basename):
    job = Job([list], [output_basename + ".snvGraphMetrics_listFiles.txt"], [['snv_graph_metrics', 'moduleVersion.R'], ['snv_graph_metrics', 'moduleVersion.tools']])

    job.command = \
"""Rscript \$R_TOOLS/snvGraphMetrics.R \\
  {list} \\
  {output_basename}""".format(
        list=list,
        output_basename=output_basename
    )

    return job

def vcf_stats(input, output, list):
    job = Job([input, list], [output], [['vcf_stats', 'moduleVersion.python'], ['vcf_stats', 'moduleVersion.tools']])

    job.command = \
"""python \$PYTHON_TOOLS/vcfStats.py \\
  -v {input} \\
  -d {dictionary} \\
  -o {output} \\
  -f {list}""".format(
        input=input,
        dictionary=config.param('vcf_stats', 'referenceSequenceDictionary', type='filepath'),
        output=output,
        list=list
    )

    return job
