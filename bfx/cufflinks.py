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
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def cuffcompare(gtf_files, output_prefix, gtf_list):
    return Job(
        gtf_files,
        [output_prefix + ".combined.gtf", output_prefix + ".TranscriptList.tsv"],
        [
            ['cuffcompare', 'module_cufflinks'],
            ['cuffcompare', 'module_mugqic_tools']
        ],
        command="""\
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
        reference_gff=config.param('cuffcompare', 'gtf', param_type='filepath'),
        reference_fasta=config.param('cuffcompare', 'genome_fasta', param_type='filepath'),
        gtf_files=" \\\n  ".join(gtf_files),
        gtf_list=gtf_list
        )
    )

def cuffdiff(sample_replicate_group_files, gtf, output_directory):

    # sample_replicate_group_files is a list of lists of replicates per sample
    # Flatten this list to set job input files
    input_files = []
    for sample_replicate_files in sample_replicate_group_files:
        input_files.extend(sample_replicate_files)

    return Job(
        input_files + [gtf],
        [os.path.join(output_directory, "isoforms.fpkm_tracking"), os.path.join(output_directory, "isoform_exp.diff")],
        [['cuffdiff', 'module_cufflinks']],
        command="""\
mkdir -p {output_directory} && \\
cuffdiff {other_options} \\
  --frag-bias-correct {genome_fasta} \\
  --library-type {library_type} \\
  --output-dir {output_directory} \\
  --num-threads {num_threads} \\
  {gtf} \\
  {input_files}""".format(
        other_options=config.param('cuffdiff', 'other_options'),
        genome_fasta=config.param('cuffdiff', 'genome_fasta', param_type='filepath'),
        library_type=config.param('cuffdiff', 'strand_info'),
        output_directory=output_directory,
        num_threads=config.param('cuffdiff', 'threads', param_type='posint'),
        gtf=gtf,
        # Join replicate bams per sample with a "," then join all sample replicate groups with a " "
        input_files=" \\\n  ".join([",".join(sample_replicate_files) for sample_replicate_files in sample_replicate_group_files])
        )
    )

def cufflinks(input_bam, output_directory, gtf=None):

    return Job(
        [input_bam],
        [os.path.join(output_directory, "transcripts.gtf"), os.path.join(output_directory, "isoforms.fpkm_tracking")],
        [['cufflinks', 'module_cufflinks']],
        command="""\
mkdir -p {output_directory} && \\
cufflinks -q {other_options}{gtf} \\
  --max-bundle-frags {max_bundle_frags} \\
  --library-type {library_type} \\
  --output-dir {output_directory} \\
  --num-threads {num_threads} \\
  {input_bam}""".format(
        other_options=config.param('cufflinks', 'other_options', required=False),
        gtf=" \\\n  --GTF-guide " + gtf if gtf else "",
        max_bundle_frags=config.param('cufflinks', 'max_bundle_frags', param_type='int'),
        library_type=config.param('cufflinks', 'strand_info'),
        output_directory=output_directory,
        num_threads=config.param('cufflinks', 'threads', param_type='posint'),
        input_bam=input_bam
        )
    )

def cuffmerge(sample_file, output_directory, gtf_file=None):

    return Job(
        [sample_file],
        [os.path.join(output_directory, "merged.gtf")],
        [['cuffmerge', 'module_cufflinks']],
        command="""\
cuffmerge {gtf} \\
  --ref-sequence {reference_sequence} \\
  -o {output_directory} \\
  --num-threads {num_threads} \\
  {sample_file}""".format(
        gtf=" \\\n  --ref-gtf " + gtf_file if gtf_file else "",
        reference_sequence=config.param('cuffmerge', 'genome_fasta', param_type='filepath', required=True),
        output_directory=output_directory,
        num_threads=config.param('cuffmerge', 'threads', param_type='posint'),
        sample_file=sample_file
        )
    )

def cuffquant(input_bam, output_directory, gtf):

    return Job(
        [input_bam,gtf],
        [os.path.join(output_directory, "abundances.cxb")],
        [['cuffquant', 'module_cufflinks']],
        command="""\
mkdir -p {output_directory} && \\
cuffquant -q {other_options} \\
  --max-bundle-frags {max_bundle_frags} \\
  --library-type {library_type} \\
  --output-dir {output_directory} \\
  --num-threads {num_threads} \\
  {gtf} \\
  {input_bam}""".format(
        other_options=config.param('cuffquant', 'other_options', required=False),
        gtf=gtf,
        max_bundle_frags=config.param('cuffquant', 'max_bundle_frags', param_type='int'),
        library_type=config.param('cuffquant', 'strand_info'),
        output_directory=output_directory,
        num_threads=config.param('cuffquant', 'threads', param_type='posint'),
        input_bam=input_bam
        )
    )

def cuffnorm(input_files, gtf, output_directory, sample_labels):

    return Job(
        input_files + [gtf],
        [os.path.join(output_directory, "isoforms.fpkm_table"), os.path.join(output_directory, "isoforms.count_table"), os.path.join(output_directory, "isoforms.attr_table"),os.path.join(output_directory, "genes.fpkm_table"), os.path.join(output_directory, "genes.count_table"), os.path.join(output_directory, "genes.attr_table")],
        [['cuffnorm', 'module_cufflinks']],
        command="""\
mkdir -p {output_directory} && \\
cuffnorm -q {other_options} \\
  --library-type {library_type} \\
  --output-dir {output_directory} \\
  --num-threads {num_threads} \\
  --labels {sample_labels} \\
  {gtf} \\
  {input_files}""".format(
        other_options=config.param('cuffnorm', 'other_options', required=False),
        gtf=gtf,
        library_type=config.param('cuffnorm', 'strand_info'),
        output_directory=output_directory,
        num_threads=config.param('cuffnorm', 'threads', param_type='posint'),
        sample_labels=sample_labels,
        input_files=" \\\n  ".join(input_files)
        )
    )
