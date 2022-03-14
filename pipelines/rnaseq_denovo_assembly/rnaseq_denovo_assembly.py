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
import glob
import logging
import os
import re
import subprocess
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import global_config_parser, _raise, SanitycheckError
from core.job import Job, concat_jobs
import utils.utils

from bfx import differential_expression
from bfx import gq_seq_utils
from bfx import rmarkdown
from bfx import samtools
from bfx import tools
from bfx import trinity
from bfx import trinotate
from bfx import blast
from bfx import exonerate
from bfx import bash_cmd as bash
from bfx import seq2fun

from pipelines import common
from pipelines.rnaseq import rnaseq


log = logging.getLogger(__name__)

class RnaSeqDeNovoAssembly(rnaseq.RnaSeqRaw):
    """
    RNA-Seq De Novo Assembly Pipeline
    =================================

    The standard Genpipes RNA-Seq De Novo Assembly pipeline now has two protocols
    to use either [Trinity](http://trinityrnaseq.sourceforge.net/)
    software suite to reconstruct transcriptomes from RNA-Seq data or Seq2Fun software suite to generate several
    informative outputs including gene abundance tables, pathway and species hit tables which
    are required for the [Networkanalyst] (https://www.networkanalyst.ca/home.xhtml). Only RNA-Seq data is used for
    both of the protocols and any reference genome or transcriptome is not required.

    The trinity De Novo Assembly pipeline, selected using the "-t trinity" parameter starts by trimming reads
    with [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
    Then reads are normalized in order to reduce memory requirement and decrease assembly runtime, using the Trinity
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

    Thus a high quality contigs assembly is created by extracting all transcripts having a functionnal annotation as defined by trinotate,
    the Top BLASTX hit and TmHMM annotations are used by default.

    Finally, different exploratory data analysis (EDA) techniques are applied to filtered isoforms expression levels.
    Main goals of expression level EDA are the detection of outliers, potential mislabeling,  to explore the homogeneity
    of biological replicates and  to appreciate the global effects of the different experimental variables.

    An HTML summary report is automatically generated by the pipeline. This report contains description of
    the sequencing experiment as well as a detailed presentation of the pipeline steps and results. Various
    Quality Control (QC) summary statistics are included in the report and additional QC analysis is accessible
    for download directly through the report. The report includes also the main references of the software and
    methods used during the analysis, together with the full list of parameters that have been passed
    to the pipeline main script.

    The Seq2Fun De Novo Assembly pipeline, selected using the "-t seq2fun" parameter directly starts with Seq2Fun
    software suit from fastq files.
    """

    def __init__(self, *args, protocol=None, **kwargs):
        self._protocol = protocol
        # Add pipeline specific arguments
        super(RnaSeqDeNovoAssembly, self).__init__(*args, **kwargs)

    @classmethod
    def argparser(cls, argparser):
        super().argparser(argparser)
        cls._argparser.add_argument("-t", "--type", help="RNAseq analysis type", dest='protocol',
                                    choices=["trinity", "seq2fun"], default="trinity")
        return cls._argparser

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
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            job = trinity.insilico_read_normalization(
                left_or_single_reads,
                right_reads,
                "fq",
                global_config_parser.param('insilico_read_normalization_readsets', 'jellyfish_memory'),
                normalization_directory,
                global_config_parser.param('insilico_read_normalization_readsets', 'cpu', required=False, param_type='int')
            )

            job.name = "insilico_read_normalization_readsets." + readset.name
            job.samples = [readset.sample]
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
                left_or_single_reads.append(os.path.join(normalization_directory, readset.name, "single.norm.fq"))
            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

        job = trinity.insilico_read_normalization(
            left_or_single_reads,
            right_reads,
            "fq",
            global_config_parser.param('insilico_read_normalization_all', 'jellyfish_memory'),
            normalization_directory_all,
            global_config_parser.param('insilico_read_normalization_all', 'cpu', required=False, param_type='int')
        )

        job.name = "insilico_read_normalization_all"
        job.samples = self.samples
        jobs.append(job)

        report_file = os.path.join("report", "RnaSeqDeNovoAssembly.insilico_read_normalization_all.md")
        normalization_stats_file = os.path.join("insilico_read_normalization", "all", "normalization.stats.tsv")
        jobs.append(
            Job(
                [normalization_stats_file, os.path.join("report", "trimReadsetTable.tsv")],
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
                name="insilico_read_normalization_all_report",
                samples=self.samples)
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

        # Trinity job
        jobs.append(
            concat_jobs([
                trinity.trinity(input_files, trinity_fasta, output_directory, reads_option),
                Job(
                    [trinity_fasta],
                    [trinity_fasta + ".zip"],
                    command="zip -j " + trinity_fasta + ".zip " + trinity_fasta
                ),
                Job(
                    [trinity_fasta],
                    [trinity_stats_prefix + ".csv", trinity_stats_prefix + ".jpg", trinity_stats_prefix + ".pdf"],
                    [['trinity', 'module_R'], ['trinity', 'module_mugqic_R_packages']],
                    command="Rscript -e 'library(gqSeqUtils); dnaFastaStats(filename = \"" + trinity_fasta + "\", type = \"trinity\", output.prefix = \"" + trinity_stats_prefix + "\")'"
                )
            ], name="trinity", samples=self.samples)
        )

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
                name="trinity_report",
                samples=self.samples)
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
        num_fasta_chunks = global_config_parser.param('exonerate_fastasplit', 'num_fasta_chunks', param_type='posint')

        return [concat_jobs([
            Job(command="rm -rf " + trinity_chunks_directory),
            Job(command="mkdir -p " + trinity_chunks_directory),
            trinity.prepare_for_blast(trinity_fasta, trinity_fasta_for_blast),
            exonerate.fastasplit(trinity_fasta_for_blast, trinity_chunks_directory, "Trinity.fa_chunk", num_fasta_chunks)
        ], name="exonerate_fastasplit.Trinity.fasta", samples=self.samples)]

    def blastx_trinity_uniprot(self):
        """
        Annotate Trinity FASTA chunks with Swiss-Prot and UniRef databases using [blastx](http://blast.ncbi.nlm.nih.gov/).
        """

        jobs = []
        trinity_chunks_directory = os.path.join("trinity_out_dir", "Trinity.fasta_chunks")
        blast_directory = "blast"
        num_fasta_chunks = global_config_parser.param('exonerate_fastasplit', 'num_fasta_chunks', param_type='posint')
        program = "blastx"
        swissprot_db = global_config_parser.param("blastx_trinity_uniprot", "swissprot_db", param_type='prefixpath')
        uniref_db = global_config_parser.param("blastx_trinity_uniprot", "uniref_db", param_type='prefixpath')
        cpu = global_config_parser.param('blastx_trinity_uniprot', 'cpu')

        # (Removed blast on uniref_db since it's too long)
        for db in [swissprot_db]:
            if not glob.glob(db + ".*phr"):
                _raise(SanitycheckError("Error: " + db + " BLAST db files do not exist!"))

            for i in range(num_fasta_chunks):
                trinity_chunk = os.path.join(trinity_chunks_directory, "Trinity.fa_chunk_{:07d}".format(i))
                query_chunk = os.path.join(blast_directory, "query_Trinity_" + os.path.basename(db) + "_chunk_{:07d}.tsv".format(i))
                blast_chunk = os.path.join(blast_directory, program + "_Trinity_" + os.path.basename(db) + "_chunk_{:07d}.tsv".format(i))
                jobs.append(
                    concat_jobs([
                        Job(command="mkdir -p " + blast_directory, removable_files=[blast_directory]),
                        Job(command="ln -s -f " + os.path.relpath(trinity_chunk, os.path.dirname(query_chunk)) + " " + query_chunk, removable_files=[blast_directory]),
                        blast.parallel_blast(trinity_chunk, query_chunk, blast_chunk, program, db, cpu),
                    ], name="blastx_trinity_uniprot." + os.path.basename(db) + ".chunk_{:07d}".format(i), samples=self.samples)
                )

        return jobs

    def blastx_trinity_uniprot_merge(self):
        """
        Merge blastx Swiss-Prot and UniRef chunks results.
        """

        jobs = []
        blast_directory = "blast"
        num_fasta_chunks = global_config_parser.param('exonerate_fastasplit', 'num_fasta_chunks', param_type='posint')
        program = "blastx"
        blast_prefix = os.path.join(blast_directory, program + "_Trinity_")
        swissprot_db = global_config_parser.param("blastx_trinity_uniprot", "swissprot_db", param_type='prefixpath')
        uniref_db = global_config_parser.param("blastx_trinity_uniprot", "uniref_db", param_type='prefixpath')

        # (Removed blast on uniref_db since it's too long)
        for db in [swissprot_db]:
            blast_chunks = [os.path.join(blast_prefix + os.path.basename(db) + "_chunk_{:07d}.tsv".format(i)) for i in range(num_fasta_chunks)]
            blast_result = os.path.join(blast_prefix + os.path.basename(db) + ".tsv")
            jobs.append(
                concat_jobs([
                    Job(
                        blast_chunks,
                        [blast_result],
                        command="cat \\\n  " + " \\\n  ".join(blast_chunks) + " \\\n  > " + blast_result
                    ),
                    Job([blast_result], [blast_result + ".zip"], command="zip -j {blast_result}.zip {blast_result}".format(blast_result=blast_result))
                ], name="blastx_trinity_" + os.path.basename(db) + "_merge", samples=self.samples)
            )

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
                name="blastx_trinity_uniprot_merge_report",
                samples=self.samples)
        )

        return jobs

    def transdecoder(self):
        """
        Identifies candidate coding regions within transcript sequences using [Transdecoder](http://transdecoder.github.io/).
        """

        trinity_fasta = os.path.join("trinity_out_dir", "Trinity.fasta")
        transdecoder_directory = os.path.join("trinotate", "transdecoder")
        transdecoder_subdirectory = os.path.join(os.path.basename(trinity_fasta) + ".transdecoder_dir")

        jobs = trinotate.transdecoder(trinity_fasta, transdecoder_directory, transdecoder_subdirectory)
        for job in jobs:
            job.samples = self.samples

        return jobs

    def hmmer(self):
        """
        Identifies protein domains using [HMMR](http://hmmer.janelia.org/).
        """

        transdecoder_directory = os.path.join("trinotate", "transdecoder")
        transdecoder_fasta = os.path.join(transdecoder_directory, "Trinity.fasta.transdecoder.pep")
        transdecoder_pfam = os.path.join(transdecoder_directory, "Trinity.fasta.transdecoder.pfam")

        jobs = trinotate.hmmer(transdecoder_directory, transdecoder_fasta, transdecoder_pfam)
        for job in jobs:
            job.samples = self.samples

        return jobs

    def rnammer_transcriptome(self):
        """
        Identify potential rRNA transcripts using [RNAmmer](http://www.cbs.dtu.dk/cgi-bin/sw_request?rnammer).
        """

        trinity_fasta = os.path.join("trinity_out_dir", "Trinity.fasta")
        rnammer_directory = os.path.join("trinotate", "rnammer")

        jobs = trinotate.rnammer_transcriptome(trinity_fasta, rnammer_directory)
        for job in jobs:
            job.samples = self.samples

        return jobs

    def blastp_transdecoder_uniprot(self):
        """
        Search Transdecoder-predicted coding regions for sequence homologies on UniProt using [blastp](http://blast.ncbi.nlm.nih.gov/).
        """

        blast_directory = os.path.join("trinotate", "blastp")
        transdecoder_fasta = os.path.join("trinotate", "transdecoder", "Trinity.fasta.transdecoder.pep")
        db = global_config_parser.param("blastp_transdecoder_uniprot", "swissprot_db", param_type='prefixpath')

        jobs = trinotate.blastp_transdecoder_uniprot(blast_directory, transdecoder_fasta, db)
        for job in jobs:
            job.samples = self.samples

        return jobs

    def signalp(self):
        """
        Predict signal peptides using [SignalP](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp).
        """

        transdecoder_fasta = os.path.join("trinotate", "transdecoder", "Trinity.fasta.transdecoder.pep")
        signalp_gff = os.path.join("trinotate", "signalp", "signalp.gff")

        jobs = trinotate.signalp(transdecoder_fasta, signalp_gff)
        for job in jobs:
            job.samples = self.samples

        return jobs

    def tmhmm(self):
        """
        Predict transmembrane regions using [TMHMM](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm).
        """

        transdecoder_fasta = os.path.join("trinotate", "transdecoder", "Trinity.fasta.transdecoder.pep")
        tmhmm_output = os.path.join("trinotate", "tmhmm", "tmhmm.out")

        jobs = trinotate.tmhmm(transdecoder_fasta, tmhmm_output)
        for job in jobs:
            job.samples = self.samples

        return jobs

    def trinotate(self):
        """
        Perform transcriptome functional annotation and analysis using [Trinotate](http://trinotate.sourceforge.net/).
        All functional annotation data is integrated into a SQLite database and a whole annotation report is created.
        """
        jobs = []

        swissprot_db = os.path.basename(global_config_parser.param("blastx_trinity_uniprot", "swissprot_db", param_type='prefixpath'))
        transdecoder_pep = os.path.join("trinotate", "transdecoder", "Trinity.fasta.transdecoder.pep")

        job = trinotate.trinotate(
            swissprot_db=swissprot_db,
            trinity_fasta=os.path.join("trinity_out_dir", "Trinity.fasta"),
            swissprot_blastx=os.path.join("blast", "blastx_Trinity_" + swissprot_db + ".tsv"),
            transdecoder_pep=transdecoder_pep,
            transdecoder_pfam=os.path.join("trinotate", "transdecoder", "Trinity.fasta.transdecoder.pfam"),
            swissprot_blastp=os.path.join("trinotate", "blastp", "blastp_" + os.path.basename(transdecoder_pep) + "_" + swissprot_db + ".tsv"),
            rnammer=os.path.join("trinotate", "rnammer", "Trinity.fasta.rnammer.gff"),
            signalp=os.path.join("trinotate", "signalp", "signalp.gff"),
            tmhmm=os.path.join("trinotate", "tmhmm", "tmhmm.out"),
            trinotate_sqlite=os.path.join("trinotate", "Trinotate.sqlite"),
            trinotate_report=os.path.join("trinotate", "trinotate_annotation_report.tsv")
        )
        job.samples = self.samples
        jobs.append(job)
        # Render Rmarkdown Report
        jobs.append(
            rmarkdown.render(
                job_input=os.path.join("trinotate", "trinotate_annotation_report.tsv"),
                job_name="trinotate_report",
                input_rmarkdown_file=os.path.join(self.report_template_dir, "RnaSeqDeNovoAssembly.trinotate.Rmd"),
                samples=self.samples,
                render_output_dir='report',
                module_section='report',
                prerun_r='report_dir="report"; source_dir="trinotate";'
            )
        )

        return jobs

    def align_and_estimate_abundance_prep_reference(self):
        """
        Index Trinity FASTA file for further abundance estimation using [Trinity align_and_estimate_abundance.pl utility](http://trinityrnaseq.sourceforge.net/analysis/abundance_estimation.html).
        """

        trinity_fasta = os.path.join("trinity_out_dir", "Trinity.fasta")
        job = trinity.align_and_estimate_abundance(trinity_fasta, output_directory="trinity_out_dir", prep_reference=True)
        job.samples = self.samples

        return [job]

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
                    _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                    "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            job = trinity.align_and_estimate_abundance(
                trinity_fasta=trinity_fasta,
                output_directory=output_directory,
                prep_reference=False,
                left_or_single_reads=left_or_single_reads,
                right_reads=right_reads,
                sample_name=sample.name
            )
            job.samples = [sample]
            jobs.append(job)

        # Generate read files and matrix of estimated abundances, send to the differential_expression directory (God bless Joel)
        output_directory = "differential_expression"

        for item in "genes", "isoforms":
            matrix = os.path.join(output_directory, item + ".counts.matrix")
            count_files = os.path.join(output_directory, item + ".counts.files")
            align_and_estimate_abundance_results = [os.path.join("align_and_estimate_abundance", sample.name, sample.name + "." + item + ".results") for sample in self.samples]
            out_prefix = os.path.join(output_directory, item)
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + os.path.join(output_directory, item)),
                Job(
                    align_and_estimate_abundance_results,
                    [count_files],
                    command="echo -e \"" + "\\n".join(align_and_estimate_abundance_results) + "\" > " + count_files,
                    samples=self.samples
                ),
                # Create isoforms and genes matrices with counts of RNA-seq fragments per feature using Trinity RSEM utility
                trinity.abundance_estimates_to_matrix(count_files, matrix, out_prefix),
                trinity.prepare_abundance_matrix_for_dge(matrix, item),
                trinity.extract_lengths_from_RSEM_output(align_and_estimate_abundance_results[0], os.path.join(output_directory, item + ".lengths.tsv"))
            ], name="align_and_estimate_abundance." + item))

        # Parse Trinotate results to obtain blast, go annotation and a filtered set of contigs
        isoforms_lengths = os.path.join(output_directory, "isoforms.lengths.tsv")
        trinotate_annotation_report = os.path.join("trinotate", "trinotate_annotation_report.tsv")
        gene_id_column = "#gene_id" if not global_config_parser.param('trinotate', 'gene_column', required=False) else global_config_parser.param('trinotate', 'gene_column', required=False)
        transcript_id_column = "transcript_id" if not global_config_parser.param('trinotate', 'transcript_column', required=False) else global_config_parser.param('trinotate', 'gene_column', required=False)
        trinotate_filters = None if not global_config_parser.param('filter_annotated_components', 'filters_trinotate', required=False) else global_config_parser.param('filter_annotated_components', 'filters_trinotate', required=False).split("\n")

        job = tools.py_parseTrinotateOutput(
            trinotate_annotation_report,
            trinotate_annotation_report + ".genes",
            trinotate_annotation_report + ".isoforms",
            gene_id_column,
            transcript_id_column,
            isoforms_lengths,
            "align_and_estimate_abundance.parse_trinotate",
            trinotate_filters
        )
        job.samples = self.samples
        jobs.append(job)

        return jobs

    def gq_seq_utils_exploratory_analysis_rnaseq_denovo(self):
        """
        Exploratory analysis using the gqSeqUtils R package.
        """

        jobs = []

        # gqSeqUtils function call
        jobs.append(concat_jobs([
            Job(command="mkdir -p exploratory", samples=self.samples),
            gq_seq_utils.exploratory_analysis_rnaseq_denovo(
                os.path.join("differential_expression", "genes.counts.matrix"),
                os.path.join("differential_expression", "genes.lengths.tsv"),
                "exploratory"
            )
        ], name="gq_seq_utils_exploratory_analysis_rnaseq_denovo"))

        # Render Rmarkdown Report
        jobs.append(
            rmarkdown.render(
                job_input=os.path.join("exploratory", "index.tsv"),
                job_name="gq_seq_utils_exploratory_analysis_rnaseq_denovo_report",
                input_rmarkdown_file=os.path.join(self.report_template_dir, "RnaSeqDeNovoAssembly.gq_seq_utils_exploratory_analysis_rnaseq.Rmd"),
                samples=self.samples,
                render_output_dir='report',
                module_section='report', # TODO: this or exploratory?
                prerun_r='report_dir="report";' # TODO: really necessary or should be hard-coded in exploratory.Rmd?
             )
        )

        return jobs


    def filter_annotated_components(self):
        """
        Filter high quality contigs based on values in trinotate annotations. Recreate a high quality contigs fasta file and run Assembly statistics using the gqSeqUtils R package.
        """

        jobs = []
        output_directory = "filtered_assembly"
        trinity_fasta = os.path.join("trinity_out_dir", "Trinity.fasta")
        trinity_filtered = os.path.join(output_directory, "Trinity.fasta")
        trinity_filtered_prefix = os.path.join(output_directory, "Trinity")
        trinity_stats_prefix = os.path.join(output_directory, "trinity_filtered.stats")
        trinotate_annotation_report_filtered = os.path.join("trinotate", "trinotate_annotation_report.tsv" + ".isoforms_filtered.tsv")

        # Use python  to extract selected headers
        jobs.append(
            concat_jobs([
                Job(command="mkdir -p " + output_directory),
                tools.py_filterAssemblyToFastaToTsv(trinity_fasta, trinotate_annotation_report_filtered, 0, trinity_filtered_prefix),
                Job(
                    [trinity_filtered],
                    [trinity_stats_prefix + ".csv", trinity_stats_prefix + ".jpg", trinity_stats_prefix + ".pdf"],
                    [['filter_annotated_components', 'module_R'], ['filter_annotated_components', 'module_mugqic_R_packages']],
                    command="""\
Rscript -e 'library(gqSeqUtils);dnaFastaStats(filename=\"{trinity_filtered}\",type=\"trinity\",output.prefix=\"{trinity_stats_prefix}\")'""".format(
    trinity_filtered=trinity_filtered,
    trinity_stats_prefix=trinity_stats_prefix
    )
                ),
                Job(
                    [trinity_filtered],
                    [trinity_filtered + ".zip"],
                    command="zip -j " + trinity_filtered + ".zip " + trinity_filtered + " " + trinity_filtered_prefix + ".tsv"
                )
            ], name="filter_annotated_components", samples=self.samples)
        )
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
    filter_string="" if not global_config_parser.param('filter_annotated_components', 'filters_trinotate', required=False) else global_config_parser.param('filter_annotated_components', 'filters_trinotate', required=False)
    ),
                name="filter_annotated_components_report",
                report_files=[report_file],
                samples=self.samples
                )
            )

        return jobs

    def gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered(self):
        """
        Exploratory analysis using the gqSeqUtils R package using a subset of filtered transcripts
        """
        # Run exploratory analysis on filtered components
        # Extract filtered components from counts file
        jobs = []
        exploratory_output_dir = os.path.join("filtered_assembly", "exploratory")
        counts_file = os.path.join("filtered_assembly", "isoforms.counts.matrix")
        trinotate_annotation_report_filtered = os.path.join("trinotate", "trinotate_annotation_report.tsv" + ".isoforms_filtered.tsv")
        trinotate_annotation_report_filtered_header = "trinotate/trinotate_annotation_report.tsv.isoforms_filtered_header.tsv"
        lengths_file = os.path.join("differential_expression", "isoforms.lengths.tsv")
        lengths_filtered_file = os.path.join("filtered_assembly", "isoforms.lengths.tsv")

        jobs.append(
            concat_jobs([
                Job(command="mkdir -p " + exploratory_output_dir),
                Job(
                    [trinotate_annotation_report_filtered],
                    [trinotate_annotation_report_filtered_header],
                    command="sed '1s/^/ \\n/' " + trinotate_annotation_report_filtered  + " > " + trinotate_annotation_report_filtered_header
                ),
                tools.py_parseMergeCsv(
                    [trinotate_annotation_report_filtered_header, os.path.join("differential_expression", "isoforms.counts.matrix")],
                    "\\\\t",
                    counts_file,
                    "\'\'",
                    left_join=True,
                    exclude="\'\'"
                ),
                tools.py_parseMergeCsv(
                    [trinotate_annotation_report_filtered_header, lengths_file],
                    "\\\\t",
                    lengths_filtered_file,
                    "\'\' transcript_id",
                    left_join=True,
                    exclude="\' \'"
                )
            ], name="filter_annotated_components_exploratory", samples=self.samples)
        )

        # gqSeqUtils function call
        jobs.append(
            concat_jobs([
                Job(command="mkdir -p " + exploratory_output_dir),
                gq_seq_utils.exploratory_analysis_rnaseq_denovo(
                    counts_file,
                    lengths_filtered_file,
                    exploratory_output_dir
                )
            ], name="gq_seq_utils_exploratory_analysis_rnaseq_denovo", samples=self.samples)
        )

        # Render Rmarkdown Report
        jobs.append(
            rmarkdown.render(
                job_input=os.path.join(exploratory_output_dir, "index.tsv"),
                job_name="gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered_report",
                input_rmarkdown_file=os.path.join(self.report_template_dir, "RnaSeqDeNovoAssembly.gq_seq_utils_exploratory_analysis_rnaseq_filtered.Rmd"),
                samples=self.samples,
                render_output_dir='report',
                module_section='report',
                prerun_r='report_dir="report/filtered_assembly"; exploratory_dir="' + exploratory_output_dir + '";'
            )
        )
        return jobs


    def differential_expression_and_goseq_rsem(self, output_directory, item, trinotate_annotation_report):
        """
        This function returns jobs related to differential gene expression analysis using [DESEQ](http://bioconductor.org/packages/release/bioc/html/DESeq.html) and [EDGER](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html).
        Merge the results of the analysis in a single csv file. Also, performs Gene Ontology analysis for RNA-Seq denovo Assembly using the Bioconductor's R package [goseq](http://www.bioconductor.org/packages/release/bioc/html/goseq.html).
        Generates GO annotations for differential genes and isoforms expression analysis, based on associated GOTERMS generated by trinotate.
        """
        jobs = []
        # Parameters from ini file
        gene_id_column = "#gene_id" if not global_config_parser.param('trinotate', 'gene_column', required=False) else global_config_parser.param('trinotate', 'gene_column', required=False)
        transcript_id_column = "transcript_id" if not global_config_parser.param('trinotate', 'transcript_column', required=False) else global_config_parser.param('trinotate', 'gene_column', required=False)
        trinotate_filters = None if not global_config_parser.param('filter_annotated_components', 'filters_trinotate', required=False) else global_config_parser.param('filter_annotated_components', 'filters_trinotate', required=False).split("\n")
        trinotate_columns_to_exclude = None if not global_config_parser.param('differential_expression', 'trinotate_columns_to_exclude', required=False) else global_config_parser.param('differential_expression', 'trinotate_columns_to_exclude', required=False)

        # mkdir
        mkdir_job = Job(
            [trinotate_annotation_report],
            [],
            command="mkdir -p " + os.path.join(output_directory, item)
        )

        # Run DGE and merge dge results with annotations
        matrix = os.path.join(output_directory, item + ".counts.matrix")

        # Perform edgeR
        edger_job = differential_expression.edger(os.path.relpath(self.design_file.name, self.output_dir), matrix + ".symbol", os.path.join(output_directory, item))
        edger_job.output_files = [os.path.join(output_directory, item, contrast.name, "edger_results.csv") for contrast in self.contrasts]

        # Perform DESeq
        deseq_job = differential_expression.deseq2(os.path.relpath(self.design_file.name, self.output_dir), matrix + ".symbol", os.path.join(output_directory, item))
        deseq_job.output_files = [os.path.join(output_directory, item, contrast.name, "dge_results.csv") for contrast in self.contrasts]

        jobs.append(
            concat_jobs([
                mkdir_job,
                edger_job,
                deseq_job,
            ], name="differential_expression.run." + item, samples=self.samples)
        )
        for contrast in self.contrasts:
            # Prepare GOseq job
            goseq_job = differential_expression.goseq(
                os.path.join(output_directory, item, contrast.name, "dge_trinotate_results.csv"),
                global_config_parser.param("differential_expression", "dge_input_columns"),
                os.path.join(output_directory, item, contrast.name, "gene_ontology_results.csv"),
                os.path.join(output_directory, item +".lengths.tsv.noheader.tsv"),
                trinotate_annotation_report + "." + item + "_go.tsv"
            )
            # Merge with annotations
            jobs.append(
                concat_jobs([
                    tools.py_parseMergeCsv(
                        [os.path.join(output_directory, item, contrast.name, "dge_results.csv"), trinotate_annotation_report + "." + item + "_blast.tsv"],
                        "\\\\t",
                        os.path.join(output_directory, item, contrast.name, "dge_trinotate_results.csv"),
                        "id " + "\"" + gene_id_column + "\"" if item == "genes" else "id " + transcript_id_column,
                        None,
                        trinotate_columns_to_exclude,
                        True,
                        "edger.p.value",
                        True
                    ),
                    # Run GOseq
                    goseq_job
                ], name="differential_expression.merge.annotations.goseq." + item + "." + contrast.name, samples=self.samples)
            )

        return jobs

    def differential_expression(self):
        """
        Performs differential gene expression analysis using [DESEQ](http://bioconductor.org/packages/release/bioc/html/DESeq.html) and [EDGER](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html).
        Merge the results of the analysis in a single csv file. Also, performs Gene Ontology analysis for RNA-Seq denovo Assembly using the Bioconductor's R package [goseq](http://www.bioconductor.org/packages/release/bioc/html/goseq.html).
        Generates GO annotations for differential genes and isoforms expression analysis, based on associated GOTERMS generated by trinotate.
        """
        output_directory = "differential_expression"
        jobs = []
        trinotate_annotation_report = os.path.join("trinotate", "trinotate_annotation_report.tsv")
        report_dir = 'report'
        input_rmarkdown_file = os.path.join(self.report_template_dir, "RnaSeqDeNovoAssembly.differential_expression_goseq.Rmd")

        # Run DGE and merge dge results with annotations
        for item in "genes", "isoforms":
            jobs.append(
                concat_jobs(
                    self.differential_expression_and_goseq_rsem(output_directory, item, trinotate_annotation_report),
                    name="differential_expression_" + item,
                    samples=self.samples)
            )

        output_files = []
        for job in jobs:
            output_files.extend([output_file for output_file in job.output_files if output_file not in output_files])

        # DGE Report
        # Render Rmarkdown Report
        jobs.append(
            rmarkdown.render(
                job_input=output_files,
                job_name="differential_expression_goseq_rnaseq_denovo_report",
                input_rmarkdown_file=input_rmarkdown_file,
                samples=self.samples,
                render_output_dir='report',
                module_section='report',
                prerun_r='design_file="' +
                         os.path.relpath(self.design_file.name, self.output_dir) +
                         '"; report_dir="' +
                         report_dir +
                         '"; source_dir="' +
                         output_directory +
                         '"; ' +
                         'top_n_results=10; contrasts=c("' +
                         '","'.join(contrast.name for contrast in self.contrasts) + '");'
            )
        )
        return jobs

    def differential_expression_filtered(self):
        """
        Differential Expression and GOSEQ analysis based on filtered transcripts and genes
        """
        output_directory = os.path.join("filtered_assembly","differential_expression")
        jobs = []
        trinotate_annotation_report = os.path.join("trinotate", "trinotate_annotation_report.tsv")
        report_dir = os.path.join("report", "filtered_assembly")
        input_rmarkdown_file = os.path.join(self.report_template_dir, "RnaSeqDeNovoAssembly.differential_expression_goseq_filtered.Rmd")

        # Filter input files
        trinotate_annotation_report_filtered = trinotate_annotation_report + ".isoforms_filtered.tsv"
        trinotate_annotation_report_filtered_header = {}
        trinotate_annotation_report_filtered_header["isoforms"] = trinotate_annotation_report + ".isoforms_filtered_header.tsv"
        trinotate_annotation_report_filtered_header["genes"] = trinotate_annotation_report + ".genes_filtered_header.tsv"
        counts_ids = {'genes':"Genes", 'isoforms':"Isoforms"}
        trinotate_filters = None if not global_config_parser.param('filter_annotated_components', 'filters_trinotate', required=False) else global_config_parser.param('filter_annotated_components', 'filters_trinotate', required=False).split("\n")
        source_directory = "differential_expression"

        # Create the files containing filtered isoforms and genes with headers
        jobs.append(
            concat_jobs([
                Job(command="mkdir -p " + output_directory),
                Job(
                    [trinotate_annotation_report_filtered],
                    [trinotate_annotation_report_filtered_header["genes"]],
                    command="cat " + trinotate_annotation_report_filtered + " | awk 'BEGIN{OFS=\"_\";FS=\"_\"}{print $1,$2}' | uniq | sed '1s/^/ \\n/' " + "  > " + trinotate_annotation_report_filtered_header["genes"],
                ),
                Job(
                    [trinotate_annotation_report_filtered],
                    [trinotate_annotation_report_filtered_header["isoforms"]],
                    command="sed '1s/^/ \\n/' " + trinotate_annotation_report_filtered  + " > " + trinotate_annotation_report_filtered_header["isoforms"]
                )
            ], name="differential_expression_filtered_get_trinotate", samples=self.samples)
        )

        # Run DGE and merge dge results with annotations
        for item in "genes", "isoforms":
            matrix = os.path.join(output_directory, item + ".counts.matrix.symbol")
            job = tools.py_parseMergeCsv(
                [trinotate_annotation_report_filtered_header[item], os.path.join(source_directory, item + ".counts.matrix.symbol")],
                "\\\\t",
                matrix,
                "\'\' " + counts_ids[item],
                left_join=True,
                exclude="\' \'"
            )
            jobs.append(
                concat_jobs([
                    job,
                    Job(
                        [os.path.join(source_directory, item +".lengths.tsv.noheader.tsv")],
                        [os.path.join(output_directory, item +".lengths.tsv.noheader.tsv")],
                        command="cp " + os.path.join(source_directory, item +".lengths.tsv.noheader.tsv") + " " + os.path.join(output_directory, item +".lengths.tsv.noheader.tsv")
                    ),
                    concat_jobs(self.differential_expression_and_goseq_rsem(output_directory, item, trinotate_annotation_report))
                ], name="differential_expression_filtered_" + item, samples=self.samples)
            )

        # Dependencies for report
        output_files = []
        for job in jobs:
            output_files.extend([output_file for output_file in job.output_files if output_file not in output_files])

        # DGE Report
        # Render Rmarkdown Report
        jobs.append(
            rmarkdown.render(
                job_input=output_files,
                job_name="differential_expression_goseq_rnaseq_denovo_filtered_report",
                input_rmarkdown_file=input_rmarkdown_file,
                samples=self.samples,
                render_output_dir='report',
                module_section='report',
                prerun_r='report_dir="' +
                         report_dir +
                         '"; source_dir="' +
                         output_directory +
                         '"; ' +
                         'top_n_results=10; contrasts=c("' +
                         '","'.join(contrast.name for contrast in self.contrasts) + '");'
            )
        )

        return jobs

    def merge_fastq(self):

        """
                this step is performed to merge fastq files if multiple readset files for one sample is present

                This step takes as input files:

                1. FASTQ files from the readset file if available
                2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """

        jobs = []
        output_directory = "merge_fastq"
        samples_to_merge = []
        merge_fastq_job = []
        #check whether a sample has more than one readsets. if so take the sample name in to a list
        for sample in self.samples:
            if len(sample.readsets) > 1:
                samples_to_merge.append(sample.name)

        # Then loop through those sample names and merge fastqs from read1 and read2
        for merge_sample in samples_to_merge:
            merge_readsets_fastq1 = []
            merge_readsets_fastq2 = []

            for readset in self.readsets:
                if merge_sample == readset.sample.name:
                    #check whether the current sample name match with an item in the sample names list
                    #if so then we know curent sample has a readset. then we need to add all the readsets in to a list

                    if readset.run_type == "PAIRED_END":
                        merge_readsets_fastq1.append(readset.fastq1)
                        merge_readsets_fastq2.append(readset.fastq2)

                        if readset.bam:
                            candidate_fastq1 = os.path.join(self.output_dir, "raw_reads", readset.sample.name,
                                                            readset.name + ".pair1.fastq.gz")
                            candidate_fastq2 = os.path.join(self.output_dir, "raw_reads", readset.sample.name,
                                                            readset.name + ".pair2.fastq.gz")
                            merge_readsets_fastq1.append(candidate_fastq1)
                            merge_readsets_fastq2.append(candidate_fastq2)

                    elif readset.run_type == "SINGLE_END":
                        merge_readsets_fastq1.append(readset.fastq1)

                        if readset.bam:
                            candidate_fastq1 = os.path.join(self.output_dir, "raw_reads", readset.sample.name,
                                                            readset.name + ".pair1.fastq.gz")
                            merge_readsets_fastq1.append(candidate_fastq1)
                    else:
                        _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                        "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))
            # if readset is paired end then two jobs will be created for the read pairs
            # if readset is single end then one job will be created for the forward read
            if(merge_readsets_fastq2):
                #job to merge double reads
                merge_fastq = concat_jobs([Job(command="mkdir -p " + os.path.join(output_directory, merge_sample)),
                                        Job(input_files=merge_readsets_fastq1,
                                            output_files=[os.path.join(output_directory, merge_sample,
                                                                      merge_sample + "_merged.pair1.fastq.gz")],
                                           command="""if (file {input1} | grep -q compressed ) ; then
ln -nsf {input1} {temp_out1}
else
gzip -c {input1} > {temp_out1} 
fi &&
if (file {input2} | grep -q compressed ) ; then
ln -nsf {input2} {temp_out2}
else 
gzip -c {input2} > {temp_out2}
fi &&
zcat  {temp_out1} {temp_out2} > {output} &&
rm {temp_out1} &&
rm {temp_out2}""".format(
                                               input1 = merge_readsets_fastq1[0],
                                               input2 = merge_readsets_fastq1[1],
                                               temp_out1 = os.path.join(output_directory, merge_sample, merge_sample + "_temp1.read1.fastq.gz"),
                                               temp_out2= os.path.join(output_directory, merge_sample,
                                                                      merge_sample + "_temp2.read1.fastq.gz"),
                                               output = os.path.join(output_directory, merge_sample, merge_sample + "_merged.pair1.fastq.gz")
                                           )),

                                           Job(input_files=merge_readsets_fastq2,
                                               output_files=[os.path.join(output_directory, merge_sample,
                                                                    merge_sample + "_merged.pair2.fastq.gz")],
                                           command="""if (file {input1} | grep -q compressed ) ; then
ln -nsf {input1} {temp_out1}
else
gzip -c {input1} > {temp_out1} 
fi &&
if (file {input2} | grep -q compressed ) ; then
ln -nsf {input2} {temp_out2}
else 
gzip -c {input2} > {temp_out2}
fi &&
zcat  {temp_out1} {temp_out2} > {output} &&
rm {temp_out1} &&
rm {temp_out2}""".format(
                                               input1 = merge_readsets_fastq2[0],
                                               input2 = merge_readsets_fastq2[1],
                                               temp_out1 = os.path.join(output_directory, merge_sample, merge_sample + "_temp1.read2.fastq.gz"),
                                               temp_out2=os.path.join(output_directory, merge_sample,
                                                                      merge_sample + "_temp2.read2.fastq.gz"),
                                               output = os.path.join(output_directory, merge_sample, merge_sample + "_merged.pair2.fastq.gz")
                                           ))
                                           ])
            else:
                #jobs to merge single reads
                merge_fastq = concat_jobs([Job(command="mkdir -p " + os.path.join(output_directory, merge_sample)),
                                        Job(input_files=merge_readsets_fastq1,
                                            output_files=[os.path.join(output_directory, merge_sample,
                                                                      merge_sample + "_merged.pair1.fastq.gz")],
                                           command="""if (file {input1} | grep -q compressed ) ; then
ln -nsf {input1} {temp_out1}
else
gzip -c {input1} > {temp_out1} 
fi &&
if (file {input2} | grep -q compressed ) ; then
ln -nsf {input2} {temp_out2}
else 
gzip -c {input2} > {temp_out2}
fi &&
zcat  {temp_out1} {temp_out2} > {output} &&
rm {temp_out1} &&
rm {temp_out2}""".format(
                                               input1 = merge_readsets_fastq1[0],
                                               input2 = merge_readsets_fastq1[1],
                                               temp_out1 = os.path.join(output_directory, merge_sample, merge_sample + "_temp1.read1.fastq.gz"),
                                               temp_out2= os.path.join(output_directory, merge_sample,
                                                                      merge_sample + "_temp2.read1.fastq.gz"),
                                               output = os.path.join(output_directory, merge_sample, merge_sample + "_merged.pair1.fastq.gz")
                                           ))
                                           ])

            merge_fastq_job.append(merge_fastq)
            merge_fastq_jobs = concat_jobs(merge_fastq_job)

        #if there is no multiple readsets for any sample, an empty job will be created.
        if len(samples_to_merge) != 0:
            job = concat_jobs([merge_fastq_jobs])
            job.samples = self.samples
            job.name = "merge_fastq"
            jobs.append(job)
        return jobs

    def seq2fun(self):
        """
        seq2fun

        This step takes as input files:

        1. FASTQ files from the readset file if available
        2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

        This step perform seq2fun analysis and generates output files including KO abundance table and KO mapped fastq files
        (https://www.seq2fun.ca/manual.xhtml#sect4) and (https://www.seq2fun.ca/manual.xhtml#sect20)

        For each contrast different folders and all the files for that particular contrast are
        generated. Therefore, only pairwise comparisons are possible
        (treatment and controls will be added according to the 1 and 2 in the design file)
        """
        jobs = []
        samples_to_merge = []
        merge_fastq_job = []
        output_directory = "seq2fun"
        # check whether a sample has more than one readsets. if so take the sample name in to a list
        #these samples will already have merged fastq and we need to use them for the seq2fun

        merge_fastq_dir = "merge_fastq"
        write_line_control_jobs = []
        write_line_treatment_jobs = []
        folder_jobs = []
        seq2fun_input_files = []
        seq2fun_output_samples = []

        #First we need to create sample table for seq2fun. Since there can be redundancy, at the end we need to remove
        #duplicate rows. If a merged fastq file is there, the path for the merged  fastq file will be added.

        for contrast in self.contrasts:

            output_file_contrast = os.path.join(output_directory, contrast.name, contrast.name + "_sample_table.temp.txt")
            folder_job = concat_jobs([Job(command="mkdir -p " + os.path.join(output_directory, contrast.name)),
            Job(command="touch " + output_file_contrast)])
            folder_jobs.append(folder_job)
            create_folder_jobs = concat_jobs(folder_jobs)

            for control in contrast.controls:

                for sample in self.samples:
                    if control.name == sample.name:

                        for readset in self.readsets:
                            if sample.name == readset.sample.name:
                                if readset.run_type == "PAIRED_END":
                                    if len(sample.readsets) > 1:
                                        #get merged fastq link
                                        candidate_fastq1 = os.path.join(merge_fastq_dir, sample.name,
                                                     sample.name + "_merged.pair1.fastq.gz")
                                        candidate_fastq2 = os.path.join(merge_fastq_dir, sample.name,
                                                                        sample.name + "_merged.pair2.fastq.gz")

                                    else:
                                        candidate_fastq1 = readset.fastq1
                                        candidate_fastq2  = readset.fastq2
                                        if readset.bam:
                                            candidate_fastq1 = os.path.join(self.output_dir, "raw_reads",
                                                                            readset.sample.name,
                                                                            readset.name + ".pair1.fastq.gz")
                                            candidate_fastq2 = os.path.join(self.output_dir, "raw_reads",
                                                                            readset.sample.name,
                                                                            readset.name + ".pair2.fastq.gz")

                                    seq2fun_input_files.append(candidate_fastq1)
                                    seq2fun_input_files.append(candidate_fastq2)
                                    seq2fun_output_samples.append(sample.name)
                                    write_line_job = Job(
                                        input_files= [candidate_fastq1, candidate_fastq2],
                                            command="""echo -e "{sample}\t{fastq1}\t{fastq2}\tcontrol" >> {file}""".format(
                                                sample=os.path.join(output_directory, contrast.name, sample.name),
                                                file=output_file_contrast,
                                                fastq1=candidate_fastq1,
                                                fastq2=candidate_fastq2
                                            )
                                        )
                                    write_line_control_jobs.append(write_line_job)
                                    sample_table_contrast_jobs = concat_jobs(write_line_control_jobs)
                                elif readset.run_type == "SINGLE_END":
                                    if len(sample.readsets) > 1:
                                        candidate_fastq1 = os.path.join(merge_fastq_dir, sample.name,
                                                                        sample.name + "merged.pair1.fastq.gz")
                                    else:
                                        candidate_fastq1 = readset.fastq1
                                        if readset.bam:
                                            candidate_fastq1 = os.path.join(self.output_dir, "raw_reads",
                                                                            readset.sample.name,
                                                                            readset.name + ".pair1.fastq.gz")

                                    seq2fun_input_files.append(candidate_fastq1)
                                    seq2fun_output_samples.append(sample.name)
                                    write_line_job = Job(
                                        input_files=[candidate_fastq1],
                                            command="""echo -e "{sample}\t{fastq1}\tcontrol" >> {file}""".format(
                                                sample=os.path.join(output_directory, contrast.name, sample.name),
                                                file=output_file_contrast,
                                                fastq1=candidate_fastq1
                                            )

                                        )
                                    write_line_control_jobs.append(write_line_job)
                                    sample_table_contrast_jobs = concat_jobs(write_line_control_jobs)
        for contrast in self.contrasts:
            #####treatments
            output_file_contrast = os.path.join(output_directory, contrast.name, contrast.name + "_sample_table.temp.txt")
            for treatment in contrast.treatments:

                for sample in self.samples:
                    if treatment.name == sample.name:

                        for readset in self.readsets:
                            if sample.name == readset.sample.name:
                                if readset.run_type == "PAIRED_END":
                                    if len(sample.readsets) > 1:

                                        candidate_fastq1 = os.path.join(merge_fastq_dir, sample.name,
                                                                        sample.name + "_merged.pair1.fastq.gz")
                                        candidate_fastq2 = os.path.join(merge_fastq_dir, sample.name,
                                                                        sample.name + "_merged.pair2.fastq.gz")

                                    else:
                                        candidate_fastq1 = readset.fastq1
                                        candidate_fastq2 = readset.fastq2
                                        if readset.bam:
                                            candidate_fastq1 = os.path.join(self.output_dir, "raw_reads",
                                                                            readset.sample.name,
                                                                            readset.name + ".pair1.fastq.gz")
                                            candidate_fastq2 = os.path.join(self.output_dir, "raw_reads",
                                                                            readset.sample.name,
                                                                            readset.name + ".pair2.fastq.gz")

                                    seq2fun_input_files.append(candidate_fastq1)
                                    seq2fun_input_files.append(candidate_fastq2)
                                    seq2fun_output_samples.append(sample.name)
                                    write_line_job = Job(
                                        input_files=[candidate_fastq1, candidate_fastq2],
                                        command="""echo -e "{sample}\t{fastq1}\t{fastq2}\ttreatment" >> {file}""".format(
                                            sample=os.path.join(output_directory, contrast.name, sample.name),
                                            file= output_file_contrast,
                                            fastq1=candidate_fastq1,
                                            fastq2=candidate_fastq2
                                        )
                                    )
                                    write_line_treatment_jobs.append(write_line_job)
                                    sample_table_treatment_jobs = concat_jobs(write_line_treatment_jobs)
                                elif readset.run_type == "SINGLE_END":
                                    if len(sample.readsets) > 1:
                                        candidate_fastq1 = os.path.join(merge_fastq_dir, sample.name,
                                                                        sample.name + "merged.pair1.fastq.gz")
                                    else:
                                        candidate_fastq1 = readset.fastq1
                                        if readset.bam:
                                            candidate_fastq1 = os.path.join(self.output_dir, "raw_reads",
                                                                            readset.sample.name,
                                                                            readset.name + ".pair1.fastq.gz")
                                    seq2fun_input_files.append(candidate_fastq1)
                                    seq2fun_output_samples.append(sample.name)
                                    write_line_job = Job(
                                        input_files=[candidate_fastq1],
                                        command="""echo -e "{sample}\t{fastq1}\ttreatment" >> {file}""".format(
                                            sample=os.path.join(output_directory, contrast.name, sample.name),
                                            file=output_file_contrast,
                                            fastq1=candidate_fastq1
                                        )

                                    )
                                    write_line_treatment_jobs.append(write_line_job)
                                    sample_table_treatment_jobs = concat_jobs(write_line_treatment_jobs)



        seq2fun_jobs = []
        ##remove duplicate lines in the sample table and execute seq2fun

        for contrast in self.contrasts:
            seq2fun_outputs =[]
            output_dir = os.path.join(output_directory, contrast.name)
            profiling = (global_config_parser.param('seq2fun', 'profiling'))
            seq2fun_outputs.append(output_dir + "/All_sample_KO_abundance_table.txt")
            if "--profiling" == profiling:
                seq2fun_outputs.append(output_dir + "/All_sample_KO_abundance_table_submit2networkanalyst.txt")
                seq2fun_outputs.append(output_dir + "/All_sample_pathway_table.txt")
                seq2fun_outputs.append(output_dir + "/All_samples.html")
                seq2fun_outputs.append(output_dir + "/All_sample_species_table.txt")


            input_file_contrast = os.path.join(output_directory, contrast.name, contrast.name + "_sample_table.temp.txt")
            output_file_contrast = os.path.join(output_directory, contrast.name, contrast.name + "_sample_table.txt")
            remove_duplicates = concat_jobs([Job(
                output_files= [output_file_contrast],
                command="""sort -u {input_file} > {output_file}""".format(
                    input_file=input_file_contrast,
                    output_file=output_file_contrast
                )

            ),Job(
                command="""rm {temp_file}""".format(
                    temp_file=input_file_contrast

                )

            ),
            seq2fun.processing(seq2fun_input_files, seq2fun_outputs, output_file_contrast, profiling)
            ])
            seq2fun_jobs.append(remove_duplicates)
            seq2fun_jobs = concat_jobs(seq2fun_jobs)



        job = concat_jobs([create_folder_jobs,sample_table_contrast_jobs, sample_table_treatment_jobs, seq2fun_jobs])
        job.samples = self.samples
        job.name = "seq2fun.processing"
        jobs.append(job)
        return jobs


    def seq2fun_count_matrix(self):
        """
        This step (only runs when running differential expression) creates common sample table that will be used in
        differential expression. The original seq2fun sample table and KO abundance table
        cannot be used in this case as differential expression needs all the samples together.
        The design file is directly used to find the contrasts.

        """

        jobs =[]
        merge_fastq_dir = "merge_fastq"
        seq2fun_input_files = []
        output_directory = "seq2fun/count_matrix"
        write_line_jobs = []
        output_sample_temp_file = os.path.join(output_directory,  "common_sample_table.temp.txt")
        output_sample_file = os.path.join(output_directory, "common_sample_table.txt")

        folder_job = concat_jobs([Job(command="mkdir -p " + os.path.join(output_directory)),
                                  Job(command="touch " + output_sample_temp_file)])

        for sample in self.samples:
                for readset in self.readsets:
                    if sample.name == readset.sample.name:
                        if readset.run_type == "PAIRED_END":
                            if len(sample.readsets) > 1:

                                candidate_fastq1 = os.path.join(merge_fastq_dir, sample.name,
                                                                sample.name + "_merged.pair1.fastq.gz")
                                candidate_fastq2 = os.path.join(merge_fastq_dir, sample.name,
                                                                sample.name + "_merged.pair2.fastq.gz")

                            else:
                                candidate_fastq1 = readset.fastq1
                                candidate_fastq2 = readset.fastq2
                                if readset.bam:
                                    candidate_fastq1 = os.path.join(self.output_dir, "raw_reads",
                                                                    readset.sample.name,
                                                                    readset.name + ".pair1.fastq.gz")
                                    candidate_fastq2 = os.path.join(self.output_dir, "raw_reads",
                                                                    readset.sample.name,
                                                                    readset.name + ".pair2.fastq.gz")

                            seq2fun_input_files.append(candidate_fastq1)
                            seq2fun_input_files.append(candidate_fastq2)
                            write_line_job = Job(
                                input_files=[candidate_fastq1, candidate_fastq2],
                                command="""echo -e "{sample}\t{fastq1}\t{fastq2}\tcontrol" >> {file}""".format(
                                    sample=os.path.join(output_directory, sample.name),
                                    file=output_sample_temp_file,
                                    fastq1=candidate_fastq1,
                                    fastq2=candidate_fastq2
                                )
                            )
                            write_line_jobs.append(write_line_job)
                            sample_table_contrast_jobs = concat_jobs(write_line_jobs)
                        elif readset.run_type == "SINGLE_END":
                            if len(sample.readsets) > 1:
                                candidate_fastq1 = os.path.join(merge_fastq_dir, sample.name,
                                                                sample.name + "merged.pair1.fastq.gz")
                            else:
                                candidate_fastq1 = readset.fastq1
                                if readset.bam:
                                    candidate_fastq1 = os.path.join(self.output_dir, "raw_reads",
                                                                    readset.sample.name,
                                                                    readset.name + ".pair1.fastq.gz")

                            seq2fun_input_files.append(candidate_fastq1)
                            write_line_job = Job(
                                input_files=[candidate_fastq1],
                                command="""echo -e "{sample}\t{fastq1}\tcontrol" >> {file}""".format(
                                    sample=os.path.join(output_directory, sample.name),
                                    file=output_sample_temp_file,
                                    fastq1=candidate_fastq1
                                )

                            )
                            write_line_jobs.append(write_line_job)
                            sample_table_contrast_jobs = concat_jobs(write_line_jobs)

        seq2fun_jobs = []
        ##remove duplicates and execute seq2fun
        seq2fun_outputs = []
        profiling = ""
        seq2fun_outputs.append(output_directory + "/All_sample_KO_abundance_table.txt")

        remove_duplicates = concat_jobs([
            Job(
            output_files=[output_sample_file],
            command="""sort -u {input_file} > {output_file}""".format(
                input_file=output_sample_temp_file,
                output_file=output_sample_file
            )
        ),
            Job(
            command="""rm {temp_file}""".format(
                temp_file=output_sample_temp_file
            )
        ),
            seq2fun.processing(seq2fun_input_files, seq2fun_outputs, output_sample_file, profiling)
        ])
        seq2fun_jobs.append(remove_duplicates)
        seq2fun_jobs = concat_jobs(seq2fun_jobs)


        job = concat_jobs(
            [folder_job, sample_table_contrast_jobs, seq2fun_jobs])
        job.samples = self.samples
        job.name = "seq2fun.count_matrix"
        jobs.append(job)
        return jobs

    def differential_expression_seq2fun(self):
        """
        Performs differential gene expression analysis using [DESEQ2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [EDGER](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html).
        Merge the results of the analysis in a single csv file.
        """
        jobs = []
        output_directory = "differential_expression/seq2fun"
        count_matrix_temp = os.path.join("seq2fun", "count_matrix", "All_sample_KO_abundance_table.txt")
        count_matrix = os.path.join("seq2fun", "count_matrix", "count_matrix.csv")
        prepare_matrix_job = Job(output_files=[count_matrix],
                                 input_files=[count_matrix_temp],
            command="""awk 'BEGIN{{FS=OFS="\t"}} NF--' {input} | awk '{{if(NR!=2){{print $0}} }}' | awk -v OFS="\t" '{{ print $1,$0}}' > {output}""".format(
                input = count_matrix_temp,
                output = count_matrix
            )
        )
        # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        if self.contrasts:
            design_file = os.path.relpath(self.design_file.name, self.output_dir)

        #check whether design file has any contrast with no replicates. If so deseq2 cannot handle it and analysis should be
        #skipped
        no_replicates = False

        #check whether the replicates are there in a contrast. If replicates are not found this step will be skipped
        for contrast in self.contrasts:
            if len(contrast.treatments) < 2:
                no_replicates = True
            elif len(contrast.controls) < 2 :
                no_replicates = True

        if no_replicates == False:
            edger_job = differential_expression.edger(design_file, count_matrix, output_directory)
            edger_job.input_files = [design_file, count_matrix]
            edger_job.output_files = [os.path.join(output_directory, contrast.name, "edger_results.csv") for contrast in self.contrasts]
            edger_job.samples = self.samples

            deseq_job = seq2fun.deseq2(design_file, count_matrix, output_directory)
            deseq_job.output_files = [os.path.join(output_directory, contrast.name, "dge_results.csv") for contrast in self.contrasts]


            deseq_job.samples = self.samples
            report_jobs = []
            for contrast in self.contrasts:
                report_job = Job(input_files = [os.path.join(output_directory, contrast.name, "dge_results.csv")],
                             output_files = [os.path.join("report", "differential_expression", "seq2fun", contrast.name, "dge_results.csv")],
                command = ("""mkdir -p {output_folder} && rm -f {output_file} && cp {input_file} {output_file}""").format(
                input_file = os.path.join(output_directory, contrast.name, "dge_results.csv"),
                output_file=os.path.join("report", "differential_expression", "seq2fun", contrast.name, "dge_results.csv"),
                output_folder = os.path.join("report", "differential_expression", "seq2fun", contrast.name)
                        ) )
                report_jobs.append(report_job)
                report_matrix_job = concat_jobs(report_jobs)
            jobs.extend(self.seq2fun_count_matrix())
            jobs.extend([concat_jobs([
                Job(command="mkdir -p " + output_directory),
                prepare_matrix_job,
                edger_job,
                deseq_job,
                report_matrix_job
            ], name="differential_expression.seq2fun")])
            return jobs
        else:
            log.info("remove all the contrasts without replicates and re-run the pipeline")
            return []

    def pathway_enrichment_seq2fun(self):

        """

        seq2fun pathway analysis using fgsea (https://bioconductor.org/packages/release/bioc/html/fgsea.html)
         and user provide universal pathway list as KEGG map ID. The differential KO expression results obtained
         from edgeR will be using as the input for the pathway enrichment analysis

        """
        jobs = []
        output_prefix = "seq2fun_ko_pathway"
        DGE_output_directory=  "differential_expression/seq2fun"
        profiling = (global_config_parser.param('seq2fun', 'profiling'))
        if "--profiling" == profiling:
            for contrast in self.contrasts:
                output_file = os.path.join(DGE_output_directory, contrast.name, "edger_results.csv")
                #html_file = os.path.join("seq2fun", contrast.name, "All_samples.html")
                output_dir = os.path.join("seq2fun_pathway", contrast.name)
                mkdir_job = Job(command="mkdir -p " + output_dir)
                pathway_job = seq2fun.ko_pathway_analysis(output_file, output_prefix,  output_dir)
                job = concat_jobs([mkdir_job, pathway_job])
                job.name = "seq2fun_pathway." + contrast.name
                job.samples = self.samples
                jobs.append(job)
        else:
            log.info("You should first run seq2fun with profiling on to conduct the pathway enrichment analysis... skipping")
        return jobs


    @property
    def step_list(self):
        return self.protocols()[self._protocol]

    def protocols(self):
        return { "trinity":
            [
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
            self.filter_annotated_components,
            self.gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered,
            self.differential_expression_filtered
            ], "seq2fun":
            [
            self.picard_sam_to_fastq,
            self.merge_fastq,
            self.seq2fun,
            self.differential_expression_seq2fun,
            self.pathway_enrichment_seq2fun
             ]
        }

def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    # Check if Genpipes must be ran inside a container
    utils.container_wrapper_argparse(__file__, argv)
    # Build help
    epilog = RnaSeqDeNovoAssembly.process_help(argv)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve', epilog=epilog)

    # populate the parser
    parser = RnaSeqDeNovoAssembly.argparser(parser)

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

    pipeline = RnaSeqDeNovoAssembly(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file,
                                    clean=clean, report=report, force=force, job_scheduler=job_scheduler, output_dir=output_dir,
                                    design_file=design_file, no_json=no_json, container=container,
                                    protocol=protocol)

    pipeline.submit_jobs()

if __name__ == '__main__':
    main()
