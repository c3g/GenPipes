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

# MUGQIC Modules
from core.config import *
from core.job import *

# Large RNA-Seq data sets, such as those exceeding 300M pairs, are best suited for in silico normalization prior to running Trinity, in order to reduce memory requirements and greatly improve upon runtimes
def insilico_read_normalization(left_or_single_reads, right_reads, sequence_type, jellyfish_memory, output_directory=None, cpu=None):

    normalization_stats_file = "normalization.stats.tsv"
    output_files = ["left.norm." + sequence_type, "right.norm." + sequence_type] if right_reads else ["single.norm." + sequence_type]
    if output_directory:
        output_files = [os.path.join(output_directory, output_file) for output_file in output_files]
        normalization_stats_file = os.path.join(output_directory, normalization_stats_file)

    output_files.append(normalization_stats_file)

    job = Job(
        left_or_single_reads + right_reads,
        output_files,
        [
            ['insilico_read_normalization', 'module_perl'],
            ['insilico_read_normalization', 'module_trinity']
        ],
        command="""\
insilico_read_normalization.pl {other_options} \\
  --seqType {sequence_type} \\
  --JM {jellyfish_memory} \\
  --max_cov {maximum_coverage} \\
  {left_or_single_reads}{right_reads}{output_directory}{cpu}""".format(
            other_options=config.param('insilico_read_normalization', 'other_options', required=False),
            sequence_type=sequence_type,
            jellyfish_memory=jellyfish_memory,
            maximum_coverage=config.param('insilico_read_normalization', 'maximum_coverage', param_type="int"),
            left_or_single_reads=" \\\n  ".join(["--left " + read for read in left_or_single_reads]) if right_reads else " \\\n  ".join(["--single " + read for read in left_or_single_reads]),
            right_reads="".join([" \\\n  --right " + read for read in right_reads]) if right_reads else "",
            output_directory=" \\\n  --output " + output_directory if output_directory else "",
            cpu=" \\\n  --CPU " + str(cpu) if cpu else ""
        ),
        removable_files=[output_directory if output_directory else None]
    )

    if output_directory:
        job = concat_jobs([
            Job(command="mkdir -p " + output_directory),
            job
        ])

    # Count normalized reads for stats
    job = concat_jobs([
        job,
        Job(
            command="""\
wc -l {output_file} | awk '{{print "# normalized {read_type} reads\t"$1 / 4}}' > {normalization_stats_file}""".format(
                output_file=output_files[0],
                read_type="paired" if right_reads else "single",
                normalization_stats_file=normalization_stats_file
            )
        )
    ])

    return job

# Create a de novo assembly from normalized readsets using [Trinity](http://trinityrnaseq.sourceforge.net/).
def trinity(input_files, trinity_fasta, output_directory, reads_option):

    return Job(
        input_files,
        [trinity_fasta],
        [
            ['trinity', 'module_perl'],
            ['trinity', 'module_java'],
            ['trinity', 'module_trinity'],
            ['trinity', 'module_bowtie'],
            ['trinity', 'module_samtools']
        ],
        command="""\
Trinity {other_options} \\
  --max_memory {max_memory} \\
  --CPU {cpu} \\
  {reads_option} \\
  --output {output_directory}""".format(
            other_options=config.param('trinity', 'other_options'),
            max_memory=config.param('trinity', 'max_memory'),
            cpu=config.param('trinity', 'cpu'),
            reads_option=reads_option,
            output_directory=output_directory
        ),
        removable_files = [output_directory]
    )

# Index Trinity FASTA file for further abundance estimation using [Trinity align_and_estimate_abundance.pl utility](http://trinityrnaseq.sourceforge.net/analysis/abundance_estimation.html).
# or run align_and_estimate_abundance perl script 

def align_and_estimate_abundance(trinity_fasta, output_directory=None, prep_reference=True, left_or_single_reads=None, right_reads=None, sample_name=None):
    # Prepare reference
    if prep_reference and left_or_single_reads is None and sample_name is None:
        job = Job(
            [trinity_fasta],
            [trinity_fasta + ".RSEM.transcripts.fa",
                trinity_fasta + ".RSEM.idx.fa"],
            [['align_and_estimate_abundance_prep_reference', 'module_perl'],
                ['align_and_estimate_abundance_prep_reference', 'module_bowtie'],
                ['align_and_estimate_abundance_prep_reference', 'module_samtools'],
                ['align_and_estimate_abundance_prep_reference', 'module_trinity']],
            command="""\
align_and_estimate_abundance.pl \\
  --transcripts {transcripts} \\
  --seqType fa \\
  --est_method RSEM \\
  --aln_method bowtie \\
  --trinity_mode \\
  --output_dir {output_directory} \\
  --prep_reference""".format(
                transcripts=trinity_fasta,
                output_directory=output_directory
            ),
            name="align_and_estimate_abundance_prep_reference")
    else :
        # Run abundance estimates
        job = Job(
            [trinity_fasta, trinity_fasta + ".RSEM.transcripts.fa", trinity_fasta + ".RSEM.idx.fa"] + left_or_single_reads + right_reads,
            [os.path.join(output_directory, sample_name + ".genes.results"),
                os.path.join(output_directory, sample_name + ".isoforms.results")],
            [['align_and_estimate_abundance_prep_reference', 'module_perl'],
                ['align_and_estimate_abundance_prep_reference', 'module_bowtie'],
                ['align_and_estimate_abundance_prep_reference', 'module_samtools'],
                ['align_and_estimate_abundance_prep_reference', 'module_trinity']],
            command="""\
align_and_estimate_abundance.pl {other_options} \\
  --transcripts {transcripts} \\
  --seqType fq \\
  --est_method RSEM \\
  --aln_method bowtie \\
  --trinity_mode \\
  --output_prefix {sample_name} \\
  --output_dir {output_directory} \\
  --thread_count {cpu} \\
  {left_or_single_reads}{right_reads}""".format(
                other_options=config.param('align_and_estimate_abundance', 'other_options'),
                transcripts=trinity_fasta,
                sample_name=sample_name,
                output_directory=output_directory,
                cpu=config.param('align_and_estimate_abundance', 'cpu', param_type='posint'),
                left_or_single_reads="--left " + ",".join(left_or_single_reads) if right_reads else "--single " + ",".join(left_or_single_reads),
                right_reads=" \\\n  --right " + ",".join(right_reads) if right_reads else ""
            ),
            name="align_and_estimate_abundance." + sample_name,
            removable_files=[
                os.path.join(output_directory, sample_name + ".bowtie.bam"),
                os.path.join(output_directory, sample_name + ".transcript.bam"),
                os.path.join(output_directory, sample_name + ".transcript.sorted.bam"),
                os.path.join(output_directory, sample_name + ".transcript.sorted.bam.bai"),
                os.path.join(output_directory, sample_name + ".bowtie.csorted.bam"),
                os.path.join(output_directory, sample_name + ".bowtie.csorted.bam.bai")
            ]
        )

    return job


def abundance_estimates_to_matrix(count_files, matrix, out_prefix):
    return Job(
        [count_files],
        [matrix],
        [
            ['align_and_estimate_abundance_prep_reference', 'module_perl'],
            ['align_and_estimate_abundance_prep_reference', 'module_trinity'],
            ['align_and_estimate_abundance_prep_reference', 'module_R']
        ],
        command="""\
abundance_estimates_to_matrix.pl \\
  --est_method RSEM \\
  --out_prefix {out_prefix} \\
  {align_and_estimate_abundance_results}""".format(
            out_prefix=out_prefix,
            align_and_estimate_abundance_results=count_files
        )
    )

def prepare_abundance_matrix_for_dge(matrix, item):
    return Job(
        [matrix],
        [matrix + ".symbol"],
        command="""\
awk -F '\\t' '{{OFS="\\t" ; print $1,$0}}' {matrix} | sed '1s/^\\t/{item}\\tSymbol/' \\
  > {matrix}.symbol""".format(
            matrix=matrix,
            item=item.title()
    ))

# Prepare the Trinity FASTA file for blast annotation (header longer then 1K characters makes the sequence to be excluded from blast)
def prepare_for_blast(trinity_fasta, trinity_fasta_for_blast):
    return Job(
        [trinity_fasta],
        [trinity_fasta_for_blast],
        command="""\
awk \'{{ print $1 }}\' {trinity_fasta}  > {trinity_fasta_for_blast}""".format(
            trinity_fasta=trinity_fasta,
            trinity_fasta_for_blast=trinity_fasta_for_blast
    ))

# Extract isoforms and genes length values from any one of sample abundance files
# edger.R requires a matrix with gene/isoform symbol as second column
def extract_lengths_from_RSEM_output(align_and_estimate_abundance_results, output):
    return Job(
        [align_and_estimate_abundance_results],
        [output, output+".noheader.tsv"],
        command="""\
cut -f 1,3,4 {align_and_estimate_abundance_results} > {output} &&  sed \'1d\' {output} > {output}.noheader.tsv""".format(
            align_and_estimate_abundance_results=align_and_estimate_abundance_results,
            output=output
        )
    )
