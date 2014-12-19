#!/usr/bin/env python

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

from bfx import gq_seq_utils
from pipelines import common

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
        return [job]

    def trinity(self):
        """
        Create a de novo assembly from normalized readsets using [Trinity](http://trinityrnaseq.sourceforge.net/).
        """

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

        trinity_job.removable_files = [output_directory]

        return [concat_jobs([
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
        ], name="trinity")]

    def exonerate_fastasplit(self):
        """
        Split the Trinity assembly FASTA into chunks for further parallel BLAST annotations.
        """

        trinity_directory = "trinity_out_dir"
        trinity_fasta = os.path.join(trinity_directory, "Trinity.fasta")
        trinity_chunks_directory = os.path.join(trinity_directory, "Trinity.fasta_chunks")
        num_fasta_chunks = config.param('exonerate_fastasplit', 'num_fasta_chunks', type='posint')

        return [concat_jobs([
            Job(command="rm -rf " + trinity_chunks_directory),
            Job(command="mkdir -p " + trinity_chunks_directory),
            Job(
                [trinity_fasta],
                # fastasplit creates FASTA chunk files numbered with 7 digits and padded with leading 0s
                [os.path.join(trinity_chunks_directory, "Trinity.fasta_chunk_{:07d}".format(i)) for i in range(num_fasta_chunks)],
                [['exonerate_fastasplit', 'module_exonerate']],
                command="fastasplit -f " + trinity_fasta + " -o " + trinity_chunks_directory + " -c " + str(num_fasta_chunks)
            )
        ], name="exonerate_fastasplit.Trinity.fasta")]

    def blastx_swissprot(self):
        """
        Annotate Trinity FASTA chunks with Swiss-Prot database using [blastx](http://blast.ncbi.nlm.nih.gov/).
        """

        jobs = []
        trinity_chunks_directory = os.path.join("trinity_out_dir", "Trinity.fasta_chunks")
        blast_directory = "blast"
        num_fasta_chunks = config.param('exonerate_fastasplit', 'num_fasta_chunks', type='posint')
        cpu = config.param('blastx_swissprot', 'cpu')
        program = "blastx"
        db = config.param("blastx_swissprot", "db", type='prefixpath')
        if not glob.glob(db + ".*phr"):
            raise Exception("Error: " + db + " BLAST db files do not exist!")

        for i in range(num_fasta_chunks):
            query_chunk = os.path.join(trinity_chunks_directory, "Trinity.fasta_chunk_{:07d}".format(i))
            blast_chunk = os.path.join(blast_directory, program + "_Trinity_" + os.path.basename(db) + "_chunk_{:07d}.tsv".format(i))
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + blast_directory),
                Job(
                    [query_chunk],
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
            ], name="blastx_swissprot.chunk_{:07d}".format(i)))

        return jobs

    def blastx_swissprot_merge(self):
        """
        Merge blastx Swiss-Prot chunks results.
        """

        blast_directory = "blast"
        num_fasta_chunks = config.param('exonerate_fastasplit', 'num_fasta_chunks', type='posint')
        program = "blastx"
        db = config.param("blastx_swissprot", "db", type='prefixpath')
        blast_chunks = [os.path.join(blast_directory, program + "_Trinity_" + os.path.basename(db) + "_chunk_{:07d}.tsv".format(i)) for i in range(num_fasta_chunks)]
        blast_result = os.path.join(blast_directory, program + "_Trinity_" + os.path.basename(db) + ".tsv")
        return [concat_jobs([
            Job(
                blast_chunks,
                [blast_result],
                command="cat \\\n  " + " \\\n  ".join(blast_chunks) + " \\\n  > " + blast_result
            ),
            Job(command="zip -j {blast_result}.zip {blast_result}".format(blast_result=blast_result))
        ], name="blastx_swissprot_merge")]

    def transdecoder(self):
        """
        Identifies candidate coding regions within transcript sequences using [Transdecoder](http://transdecoder.sourceforge.net/).
        """

        trinity_fasta = os.path.join("trinity_out_dir", "Trinity.fasta")
        transdecoder_directory = os.path.join("trinotate", "transdecoder")

        return [concat_jobs([
            Job(command="mkdir -p " + transdecoder_directory),
            Job(command="cd " + transdecoder_directory),
            Job(command="rm -f " + "Trinity.fasta.transdecoder.pfam.dat.domtbl"),
            Job(
                [trinity_fasta],
                [os.path.join(transdecoder_directory, "Trinity.fasta.transdecoder.pep"),
                 os.path.join(transdecoder_directory, "Trinity.fasta.transdecoder.pfam.dat.domtbl")],
                [['transdecoder', 'module_perl'],
                 ['transdecoder', 'module_cd_hit'],
                 ['transdecoder', 'module_hmmer'],
                 ['transdecoder', 'module_trinity']],
                command="""\
$TRINITY_HOME/trinity-plugins/transdecoder/TransDecoder {other_options} \\
  --t {transcripts} \\
  --search_pfam {pfam_db} \\
  --CPU {cpu}""".format(
                    other_options=config.param('transdecoder', 'other_options', required=False),
                    transcripts=os.path.relpath(trinity_fasta, transdecoder_directory),
                    pfam_db=config.param('transdecoder', 'pfam_db', type='filepath'),
                    cpu=config.param('transdecoder', 'cpu', type='posint')
                )
            ),
            Job(command="cd " + os.path.join("..", "..")),
        ], name="transdecoder")]

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

    def blastp_swissprot(self):
        """
        Search Transdecoder-predicted coding regions for sequence homologies using [blastp](http://blast.ncbi.nlm.nih.gov/).
        """

        blast_directory = os.path.join("trinotate", "blastp")
        cpu = config.param('blastp_swissprot', 'cpu')
        program = "blastp"
        db = config.param("blastp_swissprot", "db", type='prefixpath')
        if not glob.glob(db + ".*phr"):
            raise Exception("Error: " + db + " BLAST db files do not exist!")

        query = os.path.join("trinotate", "transdecoder", "Trinity.fasta.transdecoder.pep")
        blast_result = os.path.join(blast_directory, program + "_" + os.path.basename(query) + "_" + os.path.basename(db) + ".tsv")
        return [concat_jobs([
            Job(command="mkdir -p " + blast_directory),
            Job(
                [query],
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
        ], name="blastp_swissprot")]

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
signalp {other_options} \\
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

        db = os.path.basename(config.param("blastx_swissprot", "db", type='prefixpath'))
        trinity_fasta = os.path.join("trinity_out_dir", "Trinity.fasta")
        blastx = os.path.join("blast", "blastx_Trinity_" + db + ".tsv")
        transdecoder_pep = os.path.join("trinotate", "transdecoder", "Trinity.fasta.transdecoder.pep")
        transdecoder_pfam = os.path.join("trinotate", "transdecoder", "Trinity.fasta.transdecoder.pfam.dat.domtbl")
        blastp = os.path.join("trinotate", "blastp", "blastp_" + os.path.basename(transdecoder_pep) + "_" + db + ".tsv")
        rnammer = os.path.join("trinotate", "rnammer", "Trinity.fasta.rnammer.gff")
        signalp = os.path.join("trinotate", "signalp", "signalp.gff")
        tmhmm = os.path.join("trinotate", "tmhmm", "tmhmm.out")
        trinotate_sqlite = os.path.join("trinotate", "Trinotate.sqlite")
        trinotate_report = os.path.join("trinotate", "trinotate_annotation_report.tsv")

        return [concat_jobs([
            Job(command="mkdir -p trinotate"),
            Job(
                [trinity_fasta,
                 blastx,
                 transdecoder_pep,
                 transdecoder_pfam,
                 blastp,
                 rnammer,
                 signalp,
                 tmhmm],
                [trinotate_sqlite, trinotate_report],
                [['trinotate', 'module_perl'], ['trinotate', 'module_trinity'], ['trinotate', 'module_trinotate']],
                command="""\
cd trinotate && \\
cp $TRINOTATE_HOME/Trinotate.sqlite . && \\
$TRINITY_HOME/util/support_scripts/get_Trinity_gene_to_trans_map.pl {trinity_fasta} > {trinity_fasta}.gene_trans_map && \\
{trinotate_command} init \\
  --gene_trans_map {trinity_fasta}.gene_trans_map \\
  --transcript {trinity_fasta} \\
  --transdecoder_pep {transdecoder_pep} && \\
{trinotate_command} LOAD_blastx {blastx} && \\
{trinotate_command} LOAD_blastp {blastp} && \\
awk '{{$4 = "cds."$4; print}}' {transdecoder_pfam} > {transdecoder_pfam}.adj && \\
{trinotate_command} LOAD_pfam {transdecoder_pfam}.adj && \\
{trinotate_command} LOAD_tmhmm {tmhmm} && \\
{trinotate_command} LOAD_signalp {signalp} && \\
{trinotate_command} LOAD_rnammer {rnammer} && \\
{trinotate_command} report -E {evalue} --pfam_cutoff {pfam_cutoff} | cut -f 1-12 \\
  > {trinotate_report}""".format(
                    trinity_fasta=os.path.join("..", trinity_fasta),
                    trinotate_command="$TRINOTATE_HOME/Trinotate " + os.path.basename(trinotate_sqlite),
                    transdecoder_pep=os.path.join("..", transdecoder_pep),
                    blastx=os.path.join("..", blastx),
                    blastp=os.path.join("..", blastp),
                    transdecoder_pfam=os.path.join("..", transdecoder_pfam),
                    tmhmm=os.path.join("..", tmhmm),
                    signalp=os.path.join("..", signalp),
                    rnammer=os.path.join("..", rnammer),
                    evalue=config.param('trinotate', 'evalue'),
                    pfam_cutoff=config.param('trinotate', 'pfam_cutoff'),
                    trinotate_report=os.path.join("..", trinotate_report)
            )),
            Job(command="cd ..")], name="trinotate")]

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
                 ['align_and_estimate_abundance_prep_reference', 'module_rsem'],
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
                 ['align_and_estimate_abundance_prep_reference', 'module_rsem'],
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
                ), name="align_and_estimate_abundance." + sample.name))

        return jobs

    def differential_expression(self):
        """
        Performs differential gene expression analysis using [DESEQ](http://bioconductor.org/packages/release/bioc/html/DESeq.html) and [EDGER](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html).
        Merge the results of the analysis in a single csv file.
        """

        jobs = []
        output_directory = "differential_expression"

        for item in "genes","isoforms":
            matrix = os.path.join(output_directory, item + ".counts.matrix")
            align_and_estimate_abundance_results = [os.path.join("align_and_estimate_abundance", sample.name, sample.name + "." + item + ".results") for sample in self.samples]
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + os.path.join(output_directory, item)),
                # Create isoforms and genes matrices with counts of RNA-seq fragments per feature using Trinity RSEM utility
                Job(
                    align_and_estimate_abundance_results,
                    [matrix],
                    [['differential_expression', 'module_perl'],
                     ['differential_expression', 'module_trinity'],
                     ['differential_expression', 'module_R']],
                    command="""\
abundance_estimates_to_matrix.pl \\
  --est_method RSEM \\
  --out_prefix {out_prefix} \\
  {align_and_estimate_abundance_results}""".format(
                        out_prefix=os.path.join(output_directory, item),
                        align_and_estimate_abundance_results=" \\\n  ".join(align_and_estimate_abundance_results)
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
awk -F '\t' '{{OFS="\t" ; print $1,$0}}' {matrix} | sed '1s/^\t/{item}\tSymbol/' \\
  > {matrix}.symbol""".format(matrix=matrix, item=item.title())
                ),
                # Perform edgeR
                Job(
                    [matrix + ".symbol"],
                    [os.path.join(output_directory, item, contrast.name, "edger_results.csv") for contrast in self.contrasts],
                    [['differential_expression', 'module_R'],
                     ['differential_expression', 'module_mugqic_tools']],
                    command="""\
Rscript $R_TOOLS/edger.R \\
  -d {design} \\
  -c {matrix}.symbol \\
  -o {output_subdirectory}""".format(
                        design=os.path.relpath(self.args.design.name, self.output_dir),
                        matrix=matrix,
                        output_subdirectory=os.path.join(output_directory, item)
                    )
                ),
                # Perform DESeq
                Job(
                    [matrix + ".symbol"],
                    [os.path.join(output_directory, item, contrast.name, "deseq_results.csv") for contrast in self.contrasts] +
                    [os.path.join(output_directory, item, contrast.name, "dge_results.csv") for contrast in self.contrasts],
                    [['differential_expression', 'module_R'],
                     ['differential_expression', 'module_mugqic_tools']],
                    command="""\
Rscript $R_TOOLS/deseq.R \\
  -d {design} \\
  -c {matrix}.symbol \\
  -o {output_subdirectory}""".format(
                        design=os.path.relpath(self.args.design.name, self.output_dir),
                        matrix=matrix,
                        output_subdirectory=os.path.join(output_directory, item)
                    )
                )
            ], name="differential_expression." + item))
        return jobs

    def gq_seq_utils_report(self):
        """
        Generates the standard report. A summary html report contains the description of
        the sequencing experiment as well as a detailed presentation of the pipeline steps and results.
        Various Quality Control (QC) summary statistics are included in the report and additional QC analysis
        is accessible for download directly through the report. The report includes also the main references
        of the software and methods used during the analysis, together with the full list of parameters
        passed to the pipeline main script.
        """

        isoforms_lengths = os.path.join("differential_expression", "isoforms.lengths.tsv")
        trinotate_annotation_report = os.path.join("trinotate", "trinotate_annotation_report.tsv")
        job = concat_jobs([
            Job([isoforms_lengths, trinotate_annotation_report],
                [trinotate_annotation_report + ".genes"],
                command="""\
sed '1d' {isoforms_lengths} | perl -pe 's/^(c\d+_g\d+)/\\1\t\\1/' | sort -k1,1 -k3,3gr | \\
awk -F "\t" 'FNR == NR {{if (!isoform[$1]) {{isoform[$1] = $2; next}}}}{{OFS="\t"; if (isoform[$1] == $2 && !gene[$1]++ || $1 == "#gene_id") {{print}}}}' - {trinotate_annotation_report} | \\
perl -pe 's/^(\S+)\t\S+\t([^^\t]*)/\\1\t\\2\t\\2/' | \\
sed '1s/^#gene_id\tTop_BLASTX_hit/Gene\tSymbol/' \\
  > {trinotate_annotation_report}.genes""".format(
                    isoforms_lengths=isoforms_lengths,
                    trinotate_annotation_report=trinotate_annotation_report
                )
            ),
            Job([trinotate_annotation_report],
                [trinotate_annotation_report + ".isoforms"],
                command="""\
cut -f 2- {trinotate_annotation_report} | awk '!isoform[$1]++' | \\
perl -pe 's/^(\S+)\t([^^\t]*)/\\1\t\\2\t\\2/' | \\
sed '1s/^transcript_id\tTop_BLASTX_hit/Transcript\tSymbol/' \\
  > {trinotate_annotation_report}.isoforms""".format(
                    trinotate_annotation_report=trinotate_annotation_report
                )
            )
        ])

        for item in "genes","isoforms":
            for contrast in self.contrasts:
                dge_results = os.path.join("differential_expression", item, contrast.name, "dge_results.csv")
                dge_trinotate_results = re.sub("results\.csv$", "trinotate_results.csv", dge_results)
                job = concat_jobs([
                    job,
                    Job([dge_results, trinotate_annotation_report + "." + item],
                        [dge_trinotate_results],
                        command="""\
awk -F "\t" 'FNR == NR {{item[$1] = $0; next}} {{OFS="\t"; print $0, item[$1]}}' \\
  {trinotate_annotation_report_item} \\
  <(sed '1s/^id/{item}/' {dge_results}) \\
  > {dge_results}.tmp && \\
paste <(cut -f1,12 {dge_results}.tmp) <(cut -f3-10,13- {dge_results}.tmp) > {dge_trinotate_results} && \\
rm {dge_results}.tmp""".format(
                            trinotate_annotation_report_item=trinotate_annotation_report + "." + item,
                            item=item.title()[0:-1],
                            dge_results=dge_results,
                            dge_trinotate_results=dge_trinotate_results
                        )
                    )
                ])

        report_job = gq_seq_utils.report(
            [config_file.name for config_file in self.args.config],
            self.output_dir,
            "RNAseqDeNovo",
            self.output_dir
        )

        report_job.input_files = [
            "metrics/trimming.stats",
            "trinity_out_dir/Trinity.fasta",
            "trinity_out_dir/Trinity.stats.csv",
            "trinity_out_dir/Trinity.stats.pdf",
        ] + [os.path.join("differential_expression", item, contrast.name, "dge_trinotate_results.csv") for item in "genes","isoforms" for contrast in self.contrasts]

        return [concat_jobs([job, report_job], name="gq_seq_utils_report")]

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
            self.blastx_swissprot,
            self.blastx_swissprot_merge,
            self.transdecoder,
            self.rnammer_transcriptome,
            self.blastp_swissprot,
            self.signalp,
            self.tmhmm,
            self.trinotate,
            self.align_and_estimate_abundance_prep_reference,
            self.align_and_estimate_abundance,
            self.differential_expression,
            self.gq_seq_utils_report
        ]

if __name__ == '__main__':
    RnaSeqDeNovoAssembly()
