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
from bfx import blast

# Identifies candidate coding regions within transcript sequences using [Transdecoder](http://transdecoder.github.io/).
def transdecoder(trinity_fasta, transdecoder_directory, transdecoder_subdirectory):

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
                other_options=global_config_parser.param('transdecoder', 'other_options', required=False),
                transcripts=os.path.relpath(trinity_fasta, transdecoder_directory),
                transdecoder_subdirectory=transdecoder_subdirectory,
                swissprot_db=global_config_parser.param('transdecoder', 'swissprot_db', param_type='prefixpath'),
                pfam_db=global_config_parser.param('transdecoder', 'pfam_db', param_type='filepath'),
                cpu=global_config_parser.param('transdecoder', 'cpu', param_type='posint')
            )
        ),
        Job(command="cd " + os.path.join("..", "..")),
    ], name="transdecoder")]

# Identifies protein domains using [HMMR](http://hmmer.janelia.org/).
def hmmer(transdecoder_directory, transdecoder_fasta, transdecoder_pfam ):
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
                cpu=global_config_parser.param('hmmer', 'cpu', param_type='posint'),
                transdecoder_pfam=transdecoder_pfam,
                pfam_db=global_config_parser.param('hmmer', 'pfam_db', param_type='filepath'),
                transdecoder_fasta=transdecoder_fasta
            )
        ),
    ], name="hmmer")]

# Identify potential rRNA transcripts using [RNAmmer](http://www.cbs.dtu.dk/cgi-bin/sw_request?rnammer).
def rnammer_transcriptome(trinity_fasta, rnammer_directory):
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
                other_options=global_config_parser.param('rnammer_transcriptome', 'other_options', required=False),
                transcriptome=os.path.relpath(trinity_fasta, rnammer_directory)
            )
        ),
        Job(command="cd " + os.path.join("..", "..")),
    ], name="rnammer_transcriptome")]

# Search Transdecoder-predicted coding regions for sequence homologies on UniProt using [blastp](http://blast.ncbi.nlm.nih.gov/).
def blastp_transdecoder_uniprot(blast_directory, transdecoder_fasta, db):
    jobs = []
    cpu = global_config_parser.param('blastp_transdecoder_uniprot', 'cpu')
    program = "blastp"

    if not glob.glob(db + ".*phr"):
        raise Exception("Error: " + db + " BLAST db files do not exist!")

    query = os.path.join(blast_directory, os.path.basename(transdecoder_fasta) + "_" + os.path.basename(db) + ".tsv")
    blast_result = os.path.join(blast_directory, program + "_" + os.path.basename(query))
    jobs.append(concat_jobs([
        Job(command="mkdir -p " + blast_directory),
        Job(command="ln -s -f " + os.path.relpath(transdecoder_fasta, os.path.dirname(query)) + " " + query, removable_files=[blast_directory]),
        blast.parallel_blast(transdecoder_fasta, query, blast_result, program, db, cpu)        
    ], name="blastp_transdecoder_uniprot." + os.path.basename(db)))

    return jobs

# Predict signal peptides using [SignalP](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp).
def signalp(transdecoder_fasta, signalp_gff ):
    return [Job(
        [transdecoder_fasta],
        [signalp_gff],
        [['signalp', 'module_perl'], ['signalp', 'module_signalp']],
        command="""\
signalp -f short {other_options} \\
-T {tmp_directory} \\
-n {signalp_gff} \\
{transdecoder_fasta}""".format(
            other_options=global_config_parser.param('signalp', 'other_options', required=False),
            tmp_directory=os.path.dirname(signalp_gff),
            signalp_gff=signalp_gff,
            transdecoder_fasta=transdecoder_fasta
        ), name="signalp")]

# Predict transmembrane regions using [TMHMM](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm).
def tmhmm(transdecoder_fasta, tmhmm_output):

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


# Perform transcriptome functional annotation and analysis using [Trinotate](http://trinotate.sourceforge.net/). All functional annotation data is integrated into a SQLite database and a whole annotation report is created.
def trinotate(swissprot_db, trinity_fasta, swissprot_blastx, transdecoder_pep, transdecoder_pfam, swissprot_blastp, rnammer, signalp, tmhmm, trinotate_sqlite, trinotate_report):
    return concat_jobs([
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
                evalue=global_config_parser.param('trinotate', 'evalue'),
                pfam_cutoff=global_config_parser.param('trinotate', 'pfam_cutoff'),
                trinotate_report=trinotate_report
        ))], name="trinotate")
