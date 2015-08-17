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
import argparse
import glob
import logging
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.readset import *
from bfx import differential_expression
from bfx import gq_seq_utils
from pipelines import common
from bfx import rmarkdown
from bfx import samtools
from bfx import tools

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
        maximum_coverage=config.param('insilico_read_normalization', 'maximum_coverage', type="int"),
        left_or_single_reads=" \\\n  ".join(["--left " + read for read in left_or_single_reads]) if right_reads else " \\\n  ".join(["--single " + read for read in left_or_single_reads]),
        right_reads="".join([" \\\n  --right " + read for read in right_reads]) if right_reads else "",
        output_directory=" \\\n  --output " + output_directory if output_directory else "",
        cpu=" \\\n  --CPU " + str(cpu) if cpu else ""
        ),
        removable_files=[output_directory if output_directory else None]
    )

    if output_directory:
        job = concat_jobs([Job(command="mkdir -p " + output_directory), job])

    # Count normalized reads for stats
    job = concat_jobs([job, Job(command="""\
wc -l {output_file} | awk '{{print "# normalized {read_type} reads\t"$1 / 4}}' > {normalization_stats_file}""".format(
        output_file=output_files[0],
        read_type="paired" if right_reads else "single",
        normalization_stats_file=normalization_stats_file
    ))])

    return job


class RnaSeqDeNovoAssembly(common.Illumina):
    """
    RNA-Seq De Novo Assembly Pipeline
    =================================

    The standard MUGQIC RNA-Seq De Novo Assembly pipeline uses the [Trinity](http://trinityrnaseq.sourceforge.net/)
    software suite to reconstruct transcriptomes from RNA-Seq data without using any reference genome or transcriptome.

    First, reads are trimmed with [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
    and normalized in order to reduce memory requirement and decrease assembly runtime, using the Trinity
    normalization utility inspired by the [Diginorm](http://arxiv.org/abs/1203.4802) algorithm.

    Then, the transcriptome is assembled on normalized reads using the Trinity assembler. Trinity creates
    a Trinity.fasta file with a list of contigs representing the transcriptome isoforms. Those transcripts
    are grouped in components mostly representing genes.

    Components and transcripts are functionally annotated using the [Trinotate](http://trinotate.sourceforge.net/) suite.

    Gene abundance estimation for each sample has been performed using [RSEM](http://deweylab.biostat.wisc.edu/rsem/)
    (RNA-Seq by Expectation-Maximization). Differential gene expression analysis is performed using
    [DESeq](http://genomebiology.com/2010/11/10/R106) and [edgeR](http://bioinformatics.oxfordjournals.org/content/26/1/139/) R Bioconductor packages.

    The DESeq and edgeR methods model **count data** by a negative binomial distribution. The parameters of
    the distribution (mean and dispersion) are estimated from the data, i.e. from the read counts in the input files.
    Both methods compute a measure of read abundance, i.e. expression level (called *base mean* or
    *mean of normalized counts* in DESeq, and *concentration* in edgeR) for each gene and apply a hypothesis test
    to each gene to evaluate differential expression. In particular, both methods determine a p-value and
    a log2 fold change (in expression level) for each gene. The Log2 FC of EdgeR is reported in the differential gene
    results file, one file per design.

    The log2fold change is the logarithm (to basis 2) of the fold change condition from condition A to B
    (mutation or treatment are the most common conditions). A "fold change" between conditions A and B at a gene
    or transcript is normally computed as the ratio at gene or transcript of the base mean of scaled counts
    for condition B to the base mean of scaled counts for condition A. Counts are scaled by a size factor in
    a step called normalization (if the counts of non-differentially expressed genes in one sample are, on average,
    twice as high as in another,  the size factor for the first sample should be twice that of the other sample).
    Each column of the count table is then divided by the size factor for this column and the count values
    are brought to a common scale, making them comparable. See the [EdgeR vignette](http://www.bioconductor.org/packages/2.12/bioc/vignettes/edgeR/inst/doc/edgeR.pdf) for additional information on normalization approaches used in the pipeline.
    
    The differential gene analysis is followed by a Gene Ontology (GO) enrichment analysis.
    This analysis use the [goseq approach](http://bioconductor.org/packages/release/bioc/html/goseq.html).
    The goseq is based on the use of non-native GO terms resulting from trinotate annotations (see details in the section 5 of
    [the corresponding vignette](http://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf).

    An HTML summary report is automatically generated by the pipeline. This report contains description of
    the sequencing experiment as well as a detailed presentation of the pipeline steps and results. Various
    Quality Control (QC) summary statistics are included in the report and additional QC analysis is accessible
    for download directly through the report. The report includes also the main references of the software and
    methods used during the analysis, together with the full list of parameters that have been passed
    to the pipeline main script.
    """

    def __init__(self):
        # Add pipeline specific arguments
        self.argparser.add_argument("-d", "--design", help="design file", type=file)

        super(RnaSeqDeNovoAssembly, self).__init__()

    def insilico_read_normalization_readsets(self):
        """
        Normalize each readset, using the Trinity normalization utility.
        """

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
        """
        Merge all normalized readsets together and normalize the result, using the Trinity normalization utility.
        """

        jobs = []
        normalization_directory = "insilico_read_normalization"
        normalization_directory_all = os.path.join(normalization_directory, "all")
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
            normalization_directory_all,
            config.param('insilico_read_normalization_all', 'cpu', required=False, type='int')
        )

        job.name = "insilico_read_normalization_all"
        jobs.append(job)

        report_file = os.path.join("report", "RnaSeqDeNovoAssembly.insilico_read_normalization_all.md")
        normalization_stats_file = os.path.join("insilico_read_normalization", "all", "normalization.stats.tsv")
        jobs.append(
            Job(
                [normalization_stats_file],
                [report_file],
                [['insilico_read_normalization_all', 'module_pandoc']],
                command="""\
mkdir -p report && \\
sum_norm=`cut -f2 {normalization_stats_file}` && \\
normalization_table=`sed '1d' report/trimReadsetTable.tsv | LC_NUMERIC=en_CA awk -v sum_norm=$sum_norm '{{sum_trim+=$4}}END{{print sprintf("%\\47d", sum_trim)"|"sprintf("%\\47d", sum_norm)"|"sprintf("%.2f", sum_norm / sum_trim * 100)}}'` && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable read_type="{read_type}" \\
  --variable normalization_table="$normalization_table" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    read_type="Paired" if self.run_type == 'PAIRED_END' else "Single",
                    normalization_stats_file=normalization_stats_file,
                    report_file=report_file
                ),
                report_files=[report_file],
                name="insilico_read_normalization_all_report")
        )

        return jobs

    def trinity(self):
        """
        Create a de novo assembly from normalized readsets using [Trinity](http://trinityrnaseq.sourceforge.net/).
        """

        jobs = []

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
            ['trinity', 'module_perl'],
            ['trinity', 'module_java'],
            ['trinity', 'module_trinity'],
            ['trinity', 'module_bowtie'],
            ['trinity', 'module_samtools']
        ])

        trinity_job.command = """\
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
        )

        trinity_job.removable_files = [output_directory]

        jobs.append(concat_jobs([
            trinity_job,
            Job(
                [trinity_fasta],
                [trinity_fasta + ".zip"],
                command="zip -j " + trinity_fasta + ".zip " + trinity_fasta),
            Job(
                [trinity_fasta],
                [trinity_stats_prefix + ".csv", trinity_stats_prefix + ".jpg", trinity_stats_prefix + ".pdf"],
                [['trinity', 'module_R'], ['trinity', 'module_mugqic_R_packages']],
                command="Rscript -e 'library(gqSeqUtils); dnaFastaStats(filename = \"" + trinity_fasta + "\", type = \"trinity\", output.prefix = \"" + trinity_stats_prefix + "\")'")
        ], name="trinity"))

        report_file = os.path.join("report", "RnaSeqDeNovoAssembly.trinity.md")
        jobs.append(
            Job(
                [trinity_fasta + ".zip", trinity_stats_prefix + ".csv", trinity_stats_prefix + ".jpg", trinity_stats_prefix + ".pdf"],
                [report_file],
                [['trinity', 'module_pandoc']],
                command="""\
mkdir -p report && \\
cp {trinity_fasta}.zip {trinity_stats_prefix}.csv {trinity_stats_prefix}.jpg {trinity_stats_prefix}.pdf report/ && \\
assembly_table=`sed '1d' {trinity_stats_prefix}.csv | perl -pe 's/^"([^"]*)",/\\1\t/g' | grep -P "^(Nb. Transcripts|Nb. Components|Total Transcripts Length|Min. Transcript Length|Median Transcript Length|Mean Transcript Length|Max. Transcript Length|N50)" | LC_NUMERIC=en_CA awk -F"\t" '{{print $1"|"sprintf("%\\47d", $2)}}'` && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable assembly_table="$assembly_table" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    trinity_fasta=trinity_fasta,
                    trinity_stats_prefix=trinity_stats_prefix,
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="trinity_report")
        )

        return jobs

    def exonerate_fastasplit(self):
        """
        Split the Trinity assembly FASTA into chunks for further parallel BLAST annotations.
        """

        trinity_directory = "trinity_out_dir"
        trinity_fasta = os.path.join(trinity_directory, "Trinity.fasta")
        trinity_fasta_for_blast = os.path.join(trinity_directory, "Trinity.fa")
        trinity_chunks_directory = os.path.join(trinity_directory, "Trinity.fasta_chunks")
        num_fasta_chunks = config.param('exonerate_fastasplit', 'num_fasta_chunks', type='posint')

        return [concat_jobs([
            Job(command="rm -rf " + trinity_chunks_directory),
            Job(command="mkdir -p " + trinity_chunks_directory),
            Job([trinity_fasta],
                [trinity_fasta_for_blast],
                command="awk \'{print $1 }\' " + trinity_fasta  + " > " +  trinity_fasta_for_blast),
            Job(
                [trinity_fasta_for_blast],
                # fastasplit creates FASTA chunk files numbered with 7 digits and padded with leading 0s
                [os.path.join(trinity_chunks_directory, "Trinity.fa_chunk_{:07d}".format(i)) for i in range(num_fasta_chunks)],
                [['exonerate_fastasplit', 'module_exonerate']],
                command="fastasplit -f " + trinity_fasta_for_blast + " -o " + trinity_chunks_directory + " -c " + str(num_fasta_chunks)
            )
        ], name="exonerate_fastasplit.Trinity.fasta")]

    def blastx_trinity_uniprot(self):
        """
        Annotate Trinity FASTA chunks with Swiss-Prot and UniRef databases using [blastx](http://blast.ncbi.nlm.nih.gov/).
        """

        jobs = []
        trinity_chunks_directory = os.path.join("trinity_out_dir", "Trinity.fasta_chunks")
        blast_directory = "blast"
        num_fasta_chunks = config.param('exonerate_fastasplit', 'num_fasta_chunks', type='posint')
        cpu = config.param('blastx_trinity_uniprot', 'cpu')
        program = "blastx"
        swissprot_db = config.param("blastx_trinity_uniprot", "swissprot_db", type='prefixpath')
        uniref_db = config.param("blastx_trinity_uniprot", "uniref_db", type='prefixpath')

        # (Removed blast on uniref_db since it's too long)
        for db in [swissprot_db]:
            if not glob.glob(db + ".*phr"):
                raise Exception("Error: " + db + " BLAST db files do not exist!")

            for i in range(num_fasta_chunks):
                trinity_chunk = os.path.join(trinity_chunks_directory, "Trinity.fa_chunk_{:07d}".format(i))
                query_chunk = os.path.join(blast_directory, "query_Trinity_" + os.path.basename(db) + "_chunk_{:07d}.tsv".format(i))
                blast_chunk = os.path.join(blast_directory, program + "_Trinity_" + os.path.basename(db) + "_chunk_{:07d}.tsv".format(i))
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + blast_directory, removable_files=[blast_directory]),
                    Job(command="ln -s -f " + os.path.relpath(trinity_chunk, os.path.dirname(query_chunk)) + " " + query_chunk, removable_files=[blast_directory]),
                    Job(
                        [trinity_chunk],
                        [blast_chunk],
                        [['blast', 'module_perl'], ['blast', 'module_mugqic_tools'], ['blast', 'module_blast']],
                        command="""\
parallelBlast.pl \\
  -file {query_chunk} \\
  --OUT {blast_chunk} \\
  -n {cpu} \\
  --BLAST "{program} -db {db} -max_target_seqs 1 -outfmt '6 std stitle'" """.format(
                            query_chunk=query_chunk,
                            blast_chunk=blast_chunk,
                            cpu=cpu,
                            program=program,
                            db=db
                        )
                    )
                ], name="blastx_trinity_uniprot." + os.path.basename(db) + ".chunk_{:07d}".format(i)))

        return jobs

    def blastx_trinity_uniprot_merge(self):
        """
        Merge blastx Swiss-Prot and UniRef chunks results.
        """

        jobs = []
        blast_directory = "blast"
        num_fasta_chunks = config.param('exonerate_fastasplit', 'num_fasta_chunks', type='posint')
        program = "blastx"
        blast_prefix = os.path.join(blast_directory, program + "_Trinity_")
        swissprot_db = config.param("blastx_trinity_uniprot", "swissprot_db", type='prefixpath')
        uniref_db = config.param("blastx_trinity_uniprot", "uniref_db", type='prefixpath')

        # (Removed blast on uniref_db since it's too long)
        for db in [swissprot_db]:
            blast_chunks = [os.path.join(blast_prefix + os.path.basename(db) + "_chunk_{:07d}.tsv".format(i)) for i in range(num_fasta_chunks)]
            blast_result = os.path.join(blast_prefix + os.path.basename(db) + ".tsv")
            jobs.append(concat_jobs([
                Job(
                    blast_chunks,
                    [blast_result],
                    command="cat \\\n  " + " \\\n  ".join(blast_chunks) + " \\\n  > " + blast_result
                ),
                Job([blast_result], [blast_result + ".zip"], command="zip -j {blast_result}.zip {blast_result}".format(blast_result=blast_result))
            ], name="blastx_trinity_" + os.path.basename(db) + "_merge"))

        report_file = os.path.join("report", "RnaSeqDeNovoAssembly.blastx_trinity_uniprot_merge.md")
        jobs.append(
            Job(
                [blast_prefix + os.path.basename(swissprot_db) + ".tsv.zip"],
                [report_file],
                [['blastx_trinity_uniprot_merge', 'module_pandoc']],
                command="""\
mkdir -p report && \\
cp {blast_prefix}{blast_db}.tsv.zip  report/ && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable blast_db="{blast_db}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    blast_prefix=blast_prefix,
                    blast_db=os.path.basename(swissprot_db),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="blastx_trinity_uniprot_merge_report")
        )

        return jobs

    def transdecoder(self):
        """
        Identifies candidate coding regions within transcript sequences using [Transdecoder](http://transdecoder.github.io/).
        """

        trinity_fasta = os.path.join("trinity_out_dir", "Trinity.fasta")
        transdecoder_directory = os.path.join("trinotate", "transdecoder")
        transdecoder_subdirectory = os.path.join(os.path.basename(trinity_fasta) + ".transdecoder_dir")

        return [concat_jobs([
            Job(command="mkdir -p " + transdecoder_directory),
            Job(command="cd " + transdecoder_directory),
            Job(
                [trinity_fasta],
                [os.path.join(transdecoder_directory, os.path.basename(trinity_fasta) + ".transdecoder.pep")],
                [['transdecoder', 'module_perl'],
                 ['transdecoder', 'module_blast'],
                 ['transdecoder', 'module_hmmer'],
                 ['transdecoder', 'module_trinotate'],
                 ['transdecoder', 'module_transdecoder']],
                command="""\
TransDecoder.LongOrfs {other_options} \\
  -t {transcripts} && \\
blastp \\
  -query {transdecoder_subdirectory}/longest_orfs.pep \\
  -db {swissprot_db} \\
  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads {cpu} \\
  > {transdecoder_subdirectory}/longest_orfs.blastp.outfmt6 && \\
hmmscan --cpu {cpu} \\
  --domtblout {transdecoder_subdirectory}/longest_orfs.pfam.domtblout \\
  {pfam_db} \\
  {transdecoder_subdirectory}/longest_orfs.pep && \\
TransDecoder.Predict \\
  -t {transcripts} \\
  --retain_pfam_hits {transdecoder_subdirectory}/longest_orfs.pfam.domtblout \\
  --retain_blastp_hits {transdecoder_subdirectory}/longest_orfs.blastp.outfmt6""".format(
                    other_options=config.param('transdecoder', 'other_options', required=False),
                    transcripts=os.path.relpath(trinity_fasta, transdecoder_directory),
                    transdecoder_subdirectory=transdecoder_subdirectory,
                    swissprot_db=config.param('transdecoder', 'swissprot_db', type='prefixpath'),
                    pfam_db=config.param('transdecoder', 'pfam_db', type='filepath'),
                    cpu=config.param('transdecoder', 'cpu', type='posint')
                )
            ),
            Job(command="cd " + os.path.join("..", "..")),
        ], name="transdecoder")]

    def hmmer(self):
        """
        Identifies protein domains using [HMMR](http://hmmer.janelia.org/).
        """

        transdecoder_directory = os.path.join("trinotate", "transdecoder")
        transdecoder_fasta = os.path.join(transdecoder_directory, "Trinity.fasta.transdecoder.pep")
        transdecoder_pfam = os.path.join(transdecoder_directory, "Trinity.fasta.transdecoder.pfam")

        return [concat_jobs([
            Job(
                [transdecoder_fasta],
                [transdecoder_pfam],
                [['hmmer', 'module_hmmer']],
                command="""\
hmmscan --cpu {cpu} \\
  --domtblout {transdecoder_pfam} \\
  {pfam_db} \\
  {transdecoder_fasta}""".format(
                    cpu=config.param('hmmer', 'cpu', type='posint'),
                    transdecoder_pfam=transdecoder_pfam,
                    pfam_db=config.param('hmmer', 'pfam_db', type='filepath'),
                    transdecoder_fasta=transdecoder_fasta
                )
            ),
        ], name="hmmer")]

    def rnammer_transcriptome(self):
        """
        Identify potential rRNA transcripts using [RNAmmer](http://www.cbs.dtu.dk/cgi-bin/sw_request?rnammer).
        """
        
        trinity_fasta = os.path.join("trinity_out_dir", "Trinity.fasta")
        rnammer_directory = os.path.join("trinotate", "rnammer")

        return [concat_jobs([
            Job(command="mkdir -p " + rnammer_directory),
            Job(command="cd " + rnammer_directory),
            Job(
                [trinity_fasta],
                [os.path.join(rnammer_directory, "Trinity.fasta.rnammer.gff")],
                [['rnammer_transcriptome', 'module_perl'],
                 ['rnammer_transcriptome', 'module_hmmer'],
                 ['rnammer_transcriptome', 'module_rnammer'],
                 ['rnammer_transcriptome', 'module_trinity'],
                 ['rnammer_transcriptome', 'module_trinotate']],
                command="""\
$TRINOTATE_HOME/util/rnammer_support/RnammerTranscriptome.pl {other_options} \\
  --transcriptome {transcriptome} \\
  --path_to_rnammer `which rnammer`""".format(
                    other_options=config.param('rnammer_transcriptome', 'other_options', required=False),
                    transcriptome=os.path.relpath(trinity_fasta, rnammer_directory)
                )
            ),
            Job(command="cd " + os.path.join("..", "..")),
        ], name="rnammer_transcriptome")]

    def blastp_transdecoder_uniprot(self):
        """
        Search Transdecoder-predicted coding regions for sequence homologies on UniProt using [blastp](http://blast.ncbi.nlm.nih.gov/).
        """

        jobs = []
        blast_directory = os.path.join("trinotate", "blastp")
        cpu = config.param('blastp_transdecoder_uniprot', 'cpu')
        program = "blastp"
        transdecoder_fasta = os.path.join("trinotate", "transdecoder", "Trinity.fasta.transdecoder.pep")
        db = config.param("blastp_transdecoder_uniprot", "swissprot_db", type='prefixpath')

        if not glob.glob(db + ".*phr"):
            raise Exception("Error: " + db + " BLAST db files do not exist!")

        query = os.path.join(blast_directory, os.path.basename(transdecoder_fasta) + "_" + os.path.basename(db) + ".tsv")
        blast_result = os.path.join(blast_directory, program + "_" + os.path.basename(query))
        jobs.append(concat_jobs([
            Job(command="mkdir -p " + blast_directory),
            Job(command="ln -s -f " + os.path.relpath(transdecoder_fasta, os.path.dirname(query)) + " " + query, removable_files=[blast_directory]),
            Job(
                [transdecoder_fasta],
                [blast_result],
                [['blast', 'module_perl'], ['blast', 'module_mugqic_tools'], ['blast', 'module_blast']],
                command="""\
parallelBlast.pl \\
  -file {query} \\
  --OUT {blast_result} \\
  -n {cpu} \\
  --BLAST "{program} -db {db} -max_target_seqs 1 -outfmt '6 std stitle'" """.format(
                    query=query,
                    blast_result=blast_result,
                    cpu=cpu,
                    program=program,
                    db=db
                )
            )
        ], name="blastp_transdecoder_uniprot." + os.path.basename(db)))

        return jobs

    def signalp(self):
        """
        Predict signal peptides using [SignalP](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp).
        """

        transdecoder_fasta = os.path.join("trinotate", "transdecoder", "Trinity.fasta.transdecoder.pep")
        signalp_gff = os.path.join("trinotate", "signalp", "signalp.gff")

        return [Job(
            [transdecoder_fasta],
            [signalp_gff],
            [['signalp', 'module_perl'], ['signalp', 'module_signalp']],
            command="""\
signalp -f short {other_options} \\
  -T {tmp_directory} \\
  -n {signalp_gff} \\
  {transdecoder_fasta}""".format(
                other_options=config.param('signalp', 'other_options', required=False),
                tmp_directory=os.path.dirname(signalp_gff),
                signalp_gff=signalp_gff,
                transdecoder_fasta=transdecoder_fasta
            ), name="signalp")]

    def tmhmm(self):
        """
        Predict transmembrane regions using [TMHMM](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm).
        """

        transdecoder_fasta = os.path.join("trinotate", "transdecoder", "Trinity.fasta.transdecoder.pep")
        tmhmm_output = os.path.join("trinotate", "tmhmm", "tmhmm.out")

        return [concat_jobs([
            Job(command="mkdir -p " + os.path.dirname(tmhmm_output)),
            Job(
                [transdecoder_fasta],
                [tmhmm_output],
                [['tmhmm', 'module_perl'], ['tmhmm', 'module_tmhmm']],
                command="""\
tmhmm --short \\
  < {transdecoder_fasta} \\
  > {tmhmm_output}""".format(
                    transdecoder_fasta=transdecoder_fasta,
                    tmhmm_output=tmhmm_output
            ))], name="tmhmm")]

    def trinotate(self):
        """
        Perform transcriptome functional annotation and analysis using [Trinotate](http://trinotate.sourceforge.net/).
        All functional annotation data is integrated into a SQLite database and a whole annotation report is created.
        """

        swissprot_db = os.path.basename(config.param("blastx_trinity_uniprot", "swissprot_db", type='prefixpath'))
        trinity_fasta = os.path.join("trinity_out_dir", "Trinity.fasta")
        swissprot_blastx = os.path.join("blast", "blastx_Trinity_" + swissprot_db + ".tsv")
        transdecoder_pep = os.path.join("trinotate", "transdecoder", "Trinity.fasta.transdecoder.pep")
        transdecoder_pfam = os.path.join("trinotate", "transdecoder", "Trinity.fasta.transdecoder.pfam")
        swissprot_blastp = os.path.join("trinotate", "blastp", "blastp_" + os.path.basename(transdecoder_pep) + "_" + swissprot_db + ".tsv")
        rnammer = os.path.join("trinotate", "rnammer", "Trinity.fasta.rnammer.gff")
        signalp = os.path.join("trinotate", "signalp", "signalp.gff")
        tmhmm = os.path.join("trinotate", "tmhmm", "tmhmm.out")
        trinotate_sqlite = os.path.join("trinotate", "Trinotate.sqlite")
        trinotate_report = os.path.join("trinotate", "trinotate_annotation_report.tsv")

        return [concat_jobs([
            Job(command="mkdir -p trinotate"),
            Job(
                [trinity_fasta,
                 swissprot_blastx,
                 transdecoder_pep,
                 transdecoder_pfam,
                 swissprot_blastp,
                 rnammer,
                 signalp,
                 tmhmm],
                [trinotate_sqlite, trinotate_report],
                [['trinotate', 'module_perl'], ['trinotate', 'module_trinity'], ['trinotate', 'module_trinotate']],
                command="""\
cp $TRINOTATE_SQLITE {trinotate_sqlite} && \\
$TRINITY_HOME/util/support_scripts/get_Trinity_gene_to_trans_map.pl \\
  {trinity_fasta} \\
  > {trinity_fasta}.gene_trans_map && \\
Trinotate {trinotate_sqlite} init \\
  --gene_trans_map {trinity_fasta}.gene_trans_map \\
  --transcript {trinity_fasta} \\
  --transdecoder_pep {transdecoder_pep} && \\
Trinotate {trinotate_sqlite} LOAD_swissprot_blastx {swissprot_blastx} && \\
Trinotate {trinotate_sqlite} LOAD_swissprot_blastp {swissprot_blastp} && \\
Trinotate {trinotate_sqlite} LOAD_pfam {transdecoder_pfam} && \\
Trinotate {trinotate_sqlite} LOAD_tmhmm {tmhmm} && \\
Trinotate {trinotate_sqlite} LOAD_signalp {signalp} && \\
Trinotate {trinotate_sqlite} LOAD_rnammer {rnammer} && \\
Trinotate {trinotate_sqlite} report -E {evalue} --pfam_cutoff {pfam_cutoff} \\
  > {trinotate_report}""".format(
                    trinity_fasta=trinity_fasta,
                    trinotate_sqlite=trinotate_sqlite,
                    transdecoder_pep=transdecoder_pep,
                    swissprot_blastx=swissprot_blastx,
                    swissprot_blastp=swissprot_blastp,
                    transdecoder_pfam=transdecoder_pfam,
                    tmhmm=tmhmm,
                    signalp=signalp,
                    rnammer=rnammer,
                    evalue=config.param('trinotate', 'evalue'),
                    pfam_cutoff=config.param('trinotate', 'pfam_cutoff'),
                    trinotate_report=trinotate_report
            ))
    ], name="trinotate")]

    def align_and_estimate_abundance_prep_reference(self):
        """
        Index Trinity FASTA file for further abundance estimation using [Trinity align_and_estimate_abundance.pl utility](http://trinityrnaseq.sourceforge.net/analysis/abundance_estimation.html).
        """

        trinity_fasta = os.path.join("trinity_out_dir", "Trinity.fasta")

        return [Job(
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
  --prep_reference""".format(transcripts=trinity_fasta),
                name="align_and_estimate_abundance_prep_reference")]

    def align_and_estimate_abundance(self):
        """
        Estimate transcript abundance using [RSEM](http://deweylab.biostat.wisc.edu/rsem/) via
        [Trinity align_and_estimate_abundance.pl utility](http://trinityrnaseq.sourceforge.net/analysis/abundance_estimation.html).
        """

        jobs = []
        trinity_fasta = os.path.join("trinity_out_dir", "Trinity.fasta")

        for sample in self.samples:
            trim_directory = os.path.join("trim", sample.name)
            output_directory = os.path.join("align_and_estimate_abundance", sample.name)
            left_or_single_reads = []
            right_reads = []

            for readset in sample.readsets:
                if readset.run_type == "PAIRED_END":
                    left_or_single_reads.append(os.path.join(trim_directory, readset.name + ".trim.pair1.fastq.gz"))
                    right_reads.append(os.path.join(trim_directory, readset.name + ".trim.pair2.fastq.gz"))
                elif readset.run_type == "SINGLE_END":
                    left_or_single_reads.append(os.path.join(trim_directory, readset.name + ".trim.single.fastq.gz"))
                else:
                    raise Exception("Error: run type \"" + readset.run_type +
                    "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            jobs.append(Job(
                [trinity_fasta, trinity_fasta + ".RSEM.transcripts.fa", trinity_fasta + ".RSEM.idx.fa"] + left_or_single_reads + right_reads,
                [os.path.join(output_directory, sample.name + ".genes.results"),
                 os.path.join(output_directory, sample.name + ".isoforms.results")],
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
  --output_prefix {sample.name} \\
  --output_dir {output_directory} \\
  --thread_count {cpu} \\
  {left_or_single_reads}{right_reads}""".format(
                    other_options=config.param('align_and_estimate_abundance', 'other_options'),
                    transcripts=trinity_fasta,
                    sample=sample,
                    output_directory=output_directory,
                    cpu=config.param('align_and_estimate_abundance', 'cpu', type='posint'),
                    left_or_single_reads="--left " + ",".join(left_or_single_reads) if right_reads else "--single " + ",".join(left_or_single_reads),
                    right_reads=" \\\n  --right " + ",".join(right_reads) if right_reads else ""
                ),
                name="align_and_estimate_abundance." + sample.name,
                removable_files=[
                    os.path.join(output_directory, sample.name + ".bowtie.bam"),
                    os.path.join(output_directory, sample.name + ".transcript.bam"),
                    os.path.join(output_directory, sample.name + ".transcript.sorted.bam"),
                    os.path.join(output_directory, sample.name + ".transcript.sorted.bam.bai"),
                    os.path.join(output_directory, sample.name + ".bowtie.csorted.bam"),
                    os.path.join(output_directory, sample.name + ".bowtie.csorted.bam.bai")
                ]
            ))
        # Generate read files and matrix of estimated abundances, send to the differential_expression directory (God bless Joel)
        output_directory = "differential_expression"
        
        for item in "genes","isoforms":
            matrix = os.path.join(output_directory, item + ".counts.matrix")
            count_files = os.path.join(output_directory, item + ".counts.files")
            align_and_estimate_abundance_results = [os.path.join("align_and_estimate_abundance", sample.name, sample.name + "." + item + ".results") for sample in self.samples]
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + os.path.join(output_directory, item)),
                Job(
                    align_and_estimate_abundance_results,
                    [count_files],
                    command="echo -e \"" + "\\n".join(align_and_estimate_abundance_results) + "\" > " + count_files
                ),
                # Create isoforms and genes matrices with counts of RNA-seq fragments per feature using Trinity RSEM utility
                Job(
                    [count_files],
                    [matrix],
                    [['align_and_estimate_abundance_prep_reference', 'module_perl'],
                     ['align_and_estimate_abundance_prep_reference', 'module_trinity'],
                     ['align_and_estimate_abundance_prep_reference', 'module_R']],
                    command="""\
abundance_estimates_to_matrix.pl \\
  --est_method RSEM \\
  --out_prefix {out_prefix} \\
  {align_and_estimate_abundance_results}""".format(
                        out_prefix=os.path.join(output_directory, item),
                        align_and_estimate_abundance_results=count_files
                    )
                ),
                # Extract isoforms and genes length values from any one of sample abundance files
                Job(
                    [align_and_estimate_abundance_results[0]],
                    [os.path.join(output_directory, item + ".lengths.tsv")],
                    command="cut -f 1,3,4 " + align_and_estimate_abundance_results[0] + " \\\n  > " + os.path.join(output_directory, item + ".lengths.tsv")),
                # edger.R requires a matrix with gene/isoform symbol as second column
                Job(
                    [matrix],
                    [matrix + ".symbol"],
                    command="""\
awk -F '\\t' '{{OFS="\\t" ; print $1,$0}}' {matrix} | sed '1s/^\\t/{item}\\tSymbol/' \\
  > {matrix}.symbol""".format(matrix=matrix, item=item.title())
                )
            ], name="align_and_estimate_abundance." + item))
                
        # Parse Trinotate results to obtain blast, go annotation and a filtered set of contigs
        isoforms_lengths = os.path.join(output_directory, "isoforms.lengths.tsv")
        trinotate_annotation_report = os.path.join("trinotate", "trinotate_annotation_report.tsv")
        gene_id_column = "#gene_id" if not config.param('trinotate', 'gene_column', required=False) else config.param('trinotate', 'gene_column', required=False) 
        transcript_id_column = "transcript_id" if not config.param('trinotate', 'transcript_column', required=False) else config.param('trinotate', 'gene_column', required=False)
        trinotate_filters = None if not config.param('filter_annotated_components', 'filters_trinotate', required=False) else config.param('filter_annotated_components', 'filters_trinotate', required=False).split("\n")
                
        jobs.append(
                tools.py_parseTrinotateOutput(trinotate_annotation_report, trinotate_annotation_report + ".genes" , trinotate_annotation_report + ".isoforms",
                    gene_id_column, 
                    transcript_id_column,
                    isoforms_lengths, 
                    "align_and_estimate_abundance.parse_trinotate",
                    trinotate_filters)
        )        

        return jobs

    def differential_expression(self):
        """
        Performs differential gene expression analysis using [DESEQ](http://bioconductor.org/packages/release/bioc/html/DESeq.html) and [EDGER](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html).
        Merge the results of the analysis in a single csv file.
        """
        output_directory = "differential_expression"
        jobs = []
        trinotate_annotation_report = os.path.join("trinotate", "trinotate_annotation_report.tsv")
        gene_id_column = "#gene_id" if not config.param('trinotate', 'gene_column', required=False) else config.param('trinotate', 'gene_column', required=False) 
        transcript_id_column = "transcript_id" if not config.param('trinotate', 'transcript_column', required=False) else config.param('trinotate', 'gene_column', required=False)
        trinotate_columns_to_exclude = None if not config.param('differential_expression', 'trinotate_columns_to_exclude', required=False) else config.param('differential_expression', 'trinotate_columns_to_exclude', required=False)

        
        # Run DGE and merge dge results with annotations
        for item in "genes","isoforms":
            matrix = os.path.join(output_directory, item + ".counts.matrix")
            # Perform edgeR
            edger_job = differential_expression.edger(os.path.relpath(self.args.design.name, self.output_dir), matrix + ".symbol", os.path.join(output_directory, item))
            edger_job.output_files = [os.path.join(output_directory, item ,contrast.name, "edger_results.csv") for contrast in self.contrasts]
            # Perform DESeq             
            deseq_job = differential_expression.deseq(os.path.relpath(self.args.design.name, self.output_dir), matrix + ".symbol", os.path.join(output_directory, item))
            deseq_job.output_files = [os.path.join(output_directory, item ,contrast.name, "dge_results.csv") for contrast in self.contrasts]
            jobs.append(concat_jobs([
                        edger_job,
                        deseq_job,
                        ], name="differential_expression.run." + item))
            for contrast in self.contrasts:
                # Merge with annotations
                     jobs.append(concat_jobs([
                        tools.py_parseMergeCsv([os.path.join("differential_expression", item, contrast.name, "dge_results.csv"), trinotate_annotation_report + "." + item + "_blast.tsv"],
                                    "\\\\t",
                                    os.path.join("differential_expression", item, contrast.name, "dge_trinotate_results.csv"),
                                    "id " + "\"" + gene_id_column + "\"" if item == "genes" else "id " + transcript_id_column, 
                                    None,
                                    trinotate_columns_to_exclude,
                                    True,
                                    "edger.p.value",
                                    True
                                    ),
                    ], name="differential_expression.merge.annotations." + item + "." + contrast.name ))

        genes_matrix = os.path.join("differential_expression", "genes.counts.matrix")
        report_file = os.path.join("report", "RnaSeqDeNovoAssembly.differential_expression.md")                
        jobs.append(
            Job(
                [genes_matrix] +
                [os.path.join("differential_expression", "genes", contrast.name, "dge_trinotate_results.csv") for contrast in self.contrasts],
                [report_file],
                [['rnaseqc', 'module_pandoc']],
                name="differential_expression_report",                
                command="""\
mkdir -p report && \\
cp {design_file} report/design.tsv && \\
sed '1s/^\\t/Gene\\t/' {genes_matrix} > report/rawCountMatrix.csv && \\
pandoc \\
  {report_template_dir}/{basename_report_file} \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable design_table="`head -7 report/design.tsv | cut -f-8 | awk -F"\t" '{{OFS="\\t"; if (NR==1) {{print; gsub(/[^\\t]/, "-")}} print}}' | sed 's/\t/|/g'`" \\
  --variable raw_count_matrix_table="`head -7 report/rawCountMatrix.csv | cut -f-8 | awk -F"\t" '{{OFS="\t"; if (NR==1) {{print; gsub(/[^\t]/, "-")}} print}}' | sed 's/|/\\\\\\\\|/g' | sed 's/\t/|/g'`" \\
  --to markdown \\
  > {report_file} && \\
for contrast in {contrasts}
do
  mkdir -p report/DiffExp/$contrast/
  echo -e "\\n#### $contrast Results\\n" >> {report_file}
  cp differential_expression/genes/$contrast/dge_trinotate_results.csv report/DiffExp/$contrast/${{contrast}}_Genes_DE_results_{swissprot_db}.tsv
  echo -e "\\nTable: Differential Gene Expression Results (**partial table**; [download full table](DiffExp/$contrast/${{contrast}}_Genes_DE_results_{swissprot_db}.tsv))\\n" >> {report_file}
  head -7 report/DiffExp/$contrast/${{contrast}}_Genes_DE_results_{swissprot_db}.tsv | cut -f-8 | sed '2i ---\t---\t---\t---\t---\t---\t---\t---' | sed 's/|/\\\\|/g' | sed 's/\t/|/g' >> {report_file}
done""".format(
                    genes_matrix=genes_matrix,
                    design_file=os.path.abspath(self.args.design.name),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file,
                    contrasts=" ".join([contrast.name for contrast in self.contrasts]),
                    swissprot_db=os.path.basename(config.param("blastx_trinity_uniprot", "swissprot_db", type='prefixpath'))
                ),
                report_files=[report_file]
        ))

        return jobs
    def gq_seq_utils_exploratory_analysis_rnaseq_denovo(self):
        """
        Exploratory analysis using the gqSeqUtils R package.
        """

        jobs = []

        # gqSeqUtils function call
        jobs.append(concat_jobs([
            Job(command="mkdir -p exploratory"),
            gq_seq_utils.exploratory_analysis_rnaseq_denovo(
                os.path.join("differential_expression", "genes.counts.matrix"),
                os.path.join("differential_expression", "genes.lengths.tsv"),
                "exploratory"
            )
        ], name="gq_seq_utils_exploratory_analysis_rnaseq_denovo"))

        # Render Rmarkdown Report
        jobs.append(
            rmarkdown.render(
             job_input            = os.path.join("exploratory", "index.tsv"),
             job_name             = "gq_seq_utils_exploratory_analysis_rnaseq_denovo_report",
             input_rmarkdown_file = os.path.join(self.report_template_dir, "RnaSeqDeNovoAssembly.gq_seq_utils_exploratory_analysis_rnaseq.Rmd") ,
             render_output_dir    = 'report',
             module_section       = 'report', # TODO: this or exploratory?
             prerun_r             = 'report_dir="report";' # TODO: really necessary or should be hard-coded in exploratory.Rmd?
             )
        )
  
        return jobs    


    def differential_expression_goseq(self):
        """
        Gene Ontology analysis for RNA-Seq denovo Assembly using the Bioconductor's R package [goseq](http://www.bioconductor.org/packages/release/bioc/html/goseq.html).
        Generates GO annotations for differential genes and isoforms expression analysis, based on associated GOTERMS generated by trinotate.
        """
        
        jobs = []
        for item in ["genes","isoforms" ]:
            jobs.append(
                Job([os.path.join("differential_expression", item +".lengths.tsv")],                  
                    [os.path.join("differential_expression", item +".lengths.tsv.noheader.tsv")],
                    [],
                    command="""sed \'1d\' {item_length_file} > {item_length_file}.noheader.tsv""".format(
                        item_length_file=os.path.join("differential_expression", item +".lengths.tsv")                        
                    ), name="format." + item + ".lengths"))
            for contrast in self.contrasts:
                # goseq for differential gene /isoforms expression results
                job = differential_expression.goseq(
                    os.path.join("differential_expression", item , contrast.name, "dge_trinotate_results.csv"),
                    config.param("differential_expression_goseq", "dge_input_columns"),
                    os.path.join("differential_expression", item, contrast.name ,"gene_ontology_results.csv"),                
                    os.path.join("differential_expression", item +".lengths.tsv.noheader.tsv"),
                    os.path.join("trinotate", "trinotate_annotation_report.tsv." + item + "_go.tsv")
                )
                job.name = "differential_expression_goseq.dge." + item + "." + contrast.name
                jobs.append(job)
        
        report_file = os.path.join("report", "RnaSeqDeNovoAssembly.differential_expression_goseq.md")
        jobs.append(
            Job(
                [os.path.join("differential_expression", "genes", contrast.name, "dge_trinotate_results.csv") for contrast in self.contrasts] +  
                [os.path.join("differential_expression", "genes", contrast.name ,"gene_ontology_results.csv") for contrast in self.contrasts],
                [report_file],
                [['differential_expression_goseq', 'module_pandoc']],
                command="""\
mkdir -p report && \\
pandoc \\
  {report_template_dir}/{basename_report_file} \\
  --template {report_template_dir}/{basename_report_file} \\
  --to markdown \\
  > {report_file} && \\
for contrast in {contrasts}
do
  mkdir -p report/DiffExp/$contrast/
  echo -e "\\n#### $contrast Results\\n" >> {report_file}
  if [ `wc -l differential_expression/genes/$contrast/gene_ontology_results.csv | cut -f1 -d\ ` -gt 1 ]
  then
    cp differential_expression/genes/$contrast/gene_ontology_results.csv report/DiffExp/$contrast/${{contrast}}_Genes_GO_results.tsv
    echo -e "\\n---\\n\\nTable: GO Results of the Differentially Expressed Genes (**partial table**; [download full table](DiffExp/${{contrast}}/${{contrast}}_Genes_GO_results.tsv))\\n" >> {report_file}
    head -7 report/DiffExp/${{contrast}}/${{contrast}}_Genes_GO_results.tsv | cut -f-8 | sed '2i ---\t---\t---\t---\t---\t---\t---\t---' | sed 's/\t/|/g' >> {report_file}  
  else
    echo -e "\\nNo FDR adjusted GO enrichment was significant (p-value too high) based on the differentially expressed gene results for this design.\\n" >> {report_file}
  fi
done""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file,
                    contrasts=" ".join([contrast.name for contrast in self.contrasts])
                ),
                report_files=[report_file],
                name="differential_expression_goseq_report")
        )
        return jobs

    def filter_annotated_components(self):
        """
        Filter high quality contigs based on values in trinotate. Recreate a high quality contigs fasta file and run Assembly statistics using the gqSeqUtils R package.
        """

        jobs = []
        output_directory = "filtered_assembly"         
        trinity_fasta = os.path.join("trinity_out_dir", "Trinity.fasta")
        trinity_filtered = os.path.join(output_directory, "Trinity.fasta")
        trinity_filtered_prefix = os.path.join(output_directory, "Trinity")
        trinity_stats_prefix = os.path.join(output_directory, "trinity_filtered.stats")
        trinotate_annotation_report_filtered = os.path.join("trinotate", "trinotate_annotation_report.tsv" + ".isoforms_filtered.tsv")
        
        # Use python  to extract selected headers, needed to force input files in mkdir job, dependencies fail using pipe_jobs (uses jobs[0].input_files and jobs[-1].output_files)
        jobs.append(concat_jobs([
            Job(command="mkdir -p " + output_directory
                ), 
            tools.py_filterAssemblyToFastaToTsv(trinity_fasta , trinotate_annotation_report_filtered, 0, trinity_filtered_prefix),
            Job(
                [trinity_filtered],
                [trinity_stats_prefix + ".csv", trinity_stats_prefix + ".jpg", trinity_stats_prefix + ".pdf"],
                [['filter_annotated_components', 'module_R'], ['filter_annotated_components', 'module_mugqic_R_packages']],
                command="""\
Rscript -e 'library(gqSeqUtils); dnaFastaStats(filename = \"{trinity_filtered}\", type = \"trinity\", output.prefix = \"{trinity_stats_prefix}\")' """.format(
                    trinity_filtered=trinity_filtered,
                    trinity_stats_prefix=trinity_stats_prefix)
            ),
            Job(
                [trinity_filtered],
                [trinity_filtered + ".zip"],
                command="zip -j " + trinity_filtered + ".zip " + trinity_filtered + " " + trinity_filtered_prefix + ".tsv"
            )
        ], name="filter_annotated_components"))
        report_file = os.path.join("report", "RnaSeqDeNovoAssembly.filtered.trinity.md")
              
        jobs.append(
                Job(
                [trinity_filtered + ".zip", trinity_stats_prefix + ".csv", trinity_stats_prefix + ".jpg", trinity_stats_prefix + ".pdf"],
                [report_file],
                [['trinity', 'module_pandoc']],
                command="""\
mkdir -p report && \\
cp {trinity_filtered}.zip report/{output_directory}.zip && \\
cp {trinity_stats_prefix}.csv {trinity_stats_prefix}.jpg {trinity_stats_prefix}.pdf report/ && \\
assembly_table=`sed '1d' {trinity_stats_prefix}.csv | perl -pe 's/^"([^"]*)",/\\1\t/g' | grep -P "^(Nb. Transcripts|Nb. Components|Total Transcripts Length|Min. Transcript Length|Median Transcript Length|Mean Transcript Length|Max. Transcript Length|N50)" | LC_NUMERIC=en_CA awk -F"\t" '{{print $1"|"sprintf("%\\47d", $2)}}'` && \\
pandoc --to=markdown \\
--template {report_template_dir}/{basename_report_file} \\
--variable assembly_table="$assembly_table" \\
--variable filter_string="{filter_string}" \\
{report_template_dir}/{basename_report_file} \\
> {report_file}""".format(
                    trinity_filtered=trinity_filtered,
                    output_directory=output_directory,
                    trinity_stats_prefix=trinity_stats_prefix,
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file,
                    filter_string="" if not config.param('filter_annotated_components', 'filters_trinotate', required=False) else config.param('filter_annotated_components', 'filters_trinotate', required=False)
                    ),                
                name="filter_annotated_components_report",
                report_files=[report_file]
                )
            )
        # Run exploratory analysis on filtered components
        # Extract filtered components from counts file    
        exploratory_output_dir = "filtered_assembly/exploratory"
        counts_file = os.path.join("filtered_assembly", "isoforms.counts.matrix")
        trinotate_annotation_report_filtered_header="trinotate/trinotate_annotation_report.tsv.isoforms_filtered_header.tsv"
        lengths_file=os.path.join("differential_expression", "isoforms.lengths.tsv")
        lengths_filtered_file = os.path.join("filtered_assembly", "isoforms.lengths.tsv")
        jobs.append(concat_jobs([
                        Job(command="mkdir -p " + exploratory_output_dir),
                        Job(command="sed '1s/^/ \\n/' " + trinotate_annotation_report_filtered  + " > " + trinotate_annotation_report_filtered_header), 
                        tools.py_parseMergeCsv([ trinotate_annotation_report_filtered_header, os.path.join("differential_expression", "isoforms.counts.matrix") ],
                                    "\\\\t",
                                    counts_file,
                                    "\'\'",
                                    left_join=True,
                                    exclude="\'\'"
                                    ),
                        tools.py_parseMergeCsv([ trinotate_annotation_report_filtered_header, lengths_file ],
                                    "\\\\t",
                                    lengths_filtered_file,
                                    "\'\' transcript_id",
                                    left_join=True,
                                    exclude="\' \'"
                                    )

                    ], name="filter_annotated_components_exploratory"))

        # gqSeqUtils function call
        jobs.append(concat_jobs([
            Job(command="mkdir -p " + exploratory_output_dir),
            gq_seq_utils.exploratory_analysis_rnaseq_denovo(
                counts_file,
                lengths_filtered_file,
                exploratory_output_dir
                
            )
        ], name="gq_seq_utils_exploratory_analysis_rnaseq_denovo"))

        # Render Rmarkdown Report
        jobs.append(
            rmarkdown.render(
             job_input            = os.path.join(exploratory_output_dir, "index.tsv"),
             job_name             = "gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_report",
             input_rmarkdown_file = os.path.join(self.report_template_dir, "RnaSeqDeNovoAssembly.gq_seq_utils_exploratory_analysis_rnaseq_filtered.Rmd") ,
             render_output_dir    = 'report',
             module_section       = 'report', 
             prerun_r             = 'report_dir="report/filtered_assembly";' 
             )
        )

        return jobs    
    
    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.insilico_read_normalization_readsets,
            self.insilico_read_normalization_all,
            self.trinity,
            self.exonerate_fastasplit,
            self.blastx_trinity_uniprot,
            self.blastx_trinity_uniprot_merge,
            self.transdecoder,
            self.hmmer,
            self.rnammer_transcriptome,
            self.blastp_transdecoder_uniprot,
            self.signalp,
            self.tmhmm,
            self.trinotate,            
            self.align_and_estimate_abundance_prep_reference,
            self.align_and_estimate_abundance,
            self.gq_seq_utils_exploratory_analysis_rnaseq_denovo,
            self.differential_expression,
            self.differential_expression_goseq,
            self.filter_annotated_components
        ]

if __name__ == '__main__':
    RnaSeqDeNovoAssembly()
