#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def dna_sample_metrics(input_directory, output, experiment_type="unknown"):
    return Job(
        [input_directory],
        [output],
        [
            ['dna_sample_metrics', 'module_R'],
            ['dna_sample_metrics', 'module_mugqic_tools']
        ],
        command="""\
Rscript $R_TOOLS/DNAsampleMetrics.R \\
  {input_directory} \\
  {output} \\
  {experiment_type}""".format(
        input_directory=input_directory,
        output=output,
        experiment_type=experiment_type
    ))

def rnaseqc(sample_file, output_directory, is_single_end=False, gtf_file=None, reference=None, ribosomal_fasta=None):
    return Job(
        [sample_file],
        [os.path.join(output_directory, "index.html"), os.path.join(output_directory, "metrics.tsv")],
        [
            ['rnaseqc', 'module_java'],
            ['rnaseqc', 'module_bwa'],
            ['rnaseqc', 'module_rnaseqc']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $RNASEQC_JAR \\
  -BWArRNA {reference_ribosomal_rna_fasta} \\
  -n {number_top_transcripts} \\
  -o {output_directory} \\
  -r {reference_genome_fasta} \\
  -s {sample_file} \\
  -t {gtf_file}{other_options}{single_end}""".format(
        tmp_dir=config.param('rnaseqc', 'tmp_dir'),
        java_other_options=config.param('rnaseqc', 'java_other_options'),
        ram=config.param('rnaseqc', 'ram'),
        reference_ribosomal_rna_fasta=ribosomal_fasta if ribosomal_fasta else config.param('rnaseqc', 'ribosomal_fasta', type='filepath'),
        number_top_transcripts=config.param('rnaseqc', 'number_top_transcript', type='int'),
        output_directory=output_directory,
        reference_genome_fasta=reference if reference else config.param('rnaseqc', 'genome_fasta', type='filepath'),
        sample_file=sample_file,
        gtf_file=gtf_file if gtf_file else config.param('rnaseqc', 'gtf', type='filepath'),
        other_options=" \\\n  " + config.param('rnaseqc', 'other_options', required=False) if config.param('rnaseqc', 'other_options', required=False) else "",
        single_end=" \\\n  -singleEnd" if is_single_end else ""
        )
    )

def rpkm_saturation(count_file, gene_size_file, rpkm_directory, saturation_directory):
    return Job(
        [count_file],
        [saturation_directory + ".zip"],
        [
            ['rpkm_saturation', 'module_R'],
            ['rpkm_saturation', 'module_mugqic_tools']
        ],
        command="""\
Rscript $R_TOOLS/rpkmSaturation.R \\
  {count_file} \\
  {gene_size_file} \\
  {rpkm_directory} \\
  {saturation_directory} \\
  {threads} \\
  {other_options} && \\
zip -r {saturation_directory}.zip {saturation_directory}""".format(
        count_file=count_file,
        gene_size_file=gene_size_file,
        rpkm_directory=rpkm_directory,
        saturation_directory=saturation_directory,
        threads=config.param('rpkm_saturation', 'threads', type='posint'),
        other_options=config.param('rpkm_saturation', 'other_options', required=False)
        ),
        removable_files=[saturation_directory]
    )

def snv_graph_metrics(list, output_basename):
    return Job(
        [list],
        [output_basename + ".snvGraphMetrics_listFiles.txt"],
        [
            ['snv_graph_metrics', 'module_R'],
            ['snv_graph_metrics', 'module_mugqic_tools']
        ],
        command="""\
Rscript $R_TOOLS/snvGraphMetrics.R \\
  {list} \\
  {output_basename}""".format(
        list=list,
        output_basename=output_basename
        )
    )

def vcf_stats(input, output, list):
    return Job(
        [input],
        [output, list],
        [
            ['vcf_stats', 'module_python'],
            ['vcf_stats', 'module_mugqic_tools']
        ],
        command="""\
python $PYTHON_TOOLS/vcfStats.py \\
  -v {input} \\
  -d {dictionary} \\
  -o {output} \\
  -f {list}""".format(
        input=input,
        dictionary=config.param('vcf_stats', 'genome_dictionary', type='filepath'),
        output=output,
        list=list
        )
    )
