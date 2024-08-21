################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes. If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import os

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job, concat_jobs

# Identifies candidate coding regions within transcript sequences using [Transdecoder](http://transdecoder.github.io/).
def transdecoder_longorfs(
    trinity_fasta,
    output_directory):

    return Job(
            [trinity_fasta],
            [os.path.join(output_directory, "longest_orfs.pep")],
            [
                ['transdecoder', 'module_perl'],
                ['transdecoder', 'module_trinotate'],
                ['transdecoder', 'module_transdecoder']
            ],
            command="""\
TransDecoder.LongOrfs {other_options} \\
  -t {transcripts} \\
  -O {output_dir}""".format(
                other_options=global_conf.global_get('transdecoder', 'longorfs_options', required=False),
                transcripts=trinity_fasta,
                output_dir=output_directory
            )
        )

def transdecoder_predict(
    trinity_fasta,
    longorf_dir,
    pfam_hits,
    blastp_hits):

    prefix = os.path.basename(trinity_fasta)
    outputs = [
            prefix + ".transdecoder.bed",
            prefix + ".transdecoder.pep",
            prefix + ".transdecoder.cds",
            prefix + ".transdecoder.gff3"
            ]

    return Job(
            [trinity_fasta, pfam_hits, blastp_hits],
            outputs,
            [
                ['transdecoder', 'module_perl'],
                ['transdecoder', 'module_trinotate'],
                ['transdecoder', 'module_transdecoder']
            ],
            command="""\
TransDecoder.Predict {other_options} \\
  -t {trinity_fasta} \\
  -O {longorf_dir} \\
  --retain_pfam_hits {pfam_hits} \\
  --retain_blastp_hits {blastp_hits}""".format(
      other_options=global_conf.global_get('transdecoder', 'predict_options', required=False),
      trinity_fasta=trinity_fasta,
      longorf_dir=longorf_dir,
      pfam_hits=pfam_hits,
      blastp_hits=blastp_hits
      )
    )

# Identifies protein domains using [HMMR](http://hmmer.janelia.org/).
def hmmer(
    transdecoder_directory, 
    transdecoder_fasta,
    transdecoder_pfam):
    
    return Job(
        [transdecoder_fasta],
        [transdecoder_pfam],
        [
            ['hmmer', 'module_hmmer']
        ],
        command="""\
hmmscan --cpu {cpu} \\
  --domtblout {transdecoder_pfam} \\
  {pfam_db} \\
  {transdecoder_fasta}""".format(
                cpu=global_conf.global_get('hmmer', 'cpu', param_type='posint'),
                transdecoder_pfam=transdecoder_pfam,
                pfam_db=global_conf.global_get('hmmer', 'pfam_db', param_type='filepath'),
                transdecoder_fasta=transdecoder_fasta
            )
        )

def parse_assembly_size(
    input
    ):

    return Job(
        [input],
        [],
        [],
        command=f"""\
export genome_size=`awk -F"," '$1 ~ "Total Transcripts Length" {{print $2/1000000 * 2}}' {input}`"""
        )

def infernal_cmscan(
    input,
    clan_file,
    cm_db,
    output,
    infernal_log,
    ini_section='infernal_cmscan'
    ):

    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_infernal']
        ],
        command="""\
cmscan -Z $genome_size {options} \\
  --tblout {output} \\
  --clanin {clan_file} \\
  {cm_db} \\
  {input} \\
  > {infernal_log}""".format(
    options=global_conf.global_get(ini_section, 'options'),
    output=output,
    clan_file=clan_file,
    cm_db=cm_db,
    input=input,
    infernal_log=infernal_log
    )
  )

def infernal_merge(
    input,
    output
    ):
    return Job(
        input,
        [output],
        [],
        command="""
head -n2 {input1} > {output}
cat {input_chunks} | grep -v '^#' >> {output}""".format(
    input1=input[0],
    input_chunks=" \\\n".join(input),
    output=output
    )
  )

def signalp(
        transdecoder_fasta,
        signalp_gff):

    return Job(
        [transdecoder_fasta],
        [signalp_gff],
        [
            ['signalp', 'module_perl'], 
            ['signalp', 'module_signalp']
        ],
        command="""\
signalp -f short {other_options} \\
  -T {tmp_directory} \\
  -n {signalp_gff} \\
  {transdecoder_fasta}""".format(
            other_options=global_conf.global_get('signalp', 'other_options', required=False),
            tmp_directory=os.path.dirname(signalp_gff),
            signalp_gff=signalp_gff,
            transdecoder_fasta=transdecoder_fasta
        )
     )

def signalp6(
        transdecoder_fasta,
        output_directory):

    output = os.path.join(output_directory, 'output.gff3')

    return Job(
        [transdecoder_fasta],
        [output],
        [
            ['signalp', 'module_python'], 
            ['signalp', 'module_signalp']
        ],
        command="""\
signalp6 {other_options} \\
  --output_dir {output_directory} \\
  --fastafile {transdecoder_fasta}""".format(
            other_options=global_conf.global_get('signalp', 'other_options', required=False),
            output_directory=output_directory,
            transdecoder_fasta=transdecoder_fasta
        )
     )

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
def trinotate(
    trinity_fasta,
    swissprot_blastx,
    transdecoder_pep,
    transdecoder_pfam,
    swissprot_blastp,
    infernal,
    signalp,
    tmhmm,
    trinotate_sqlite,
    trinotate_report
    ):
    return Job(
        [
            trinity_fasta,
            swissprot_blastx,
            transdecoder_pep,
            transdecoder_pfam,
            swissprot_blastp,
            infernal,
            signalp,
            tmhmm
        ],
        [
            trinotate_sqlite,
            trinotate_report
        ],
        [
            ['trinotate', 'module_perl'],
            ['trinotate', 'module_trinity'],
            ['trinotate', 'module_trinotate']
        ],
        command="""\
cp $TRINOTATE_SQLITE {trinotate_sqlite} && \\
Trinotate --db {trinotate_sqlite} --init \\
  --gene_trans_map {trinity_fasta}.gene_trans_map \\
  --transcript_fasta {trinity_fasta} \\
  --transdecoder_pep {transdecoder_pep} && \\
Trinotate --db {trinotate_sqlite} --LOAD_swissprot_blastx {swissprot_blastx} && \\
Trinotate --db {trinotate_sqlite} --LOAD_swissprot_blastp {swissprot_blastp} && \\
Trinotate --db {trinotate_sqlite} --LOAD_pfam {transdecoder_pfam} && \\
Trinotate --db {trinotate_sqlite} --LOAD_tmhmmv2 {tmhmm} && \\
Trinotate --db {trinotate_sqlite} --LOAD_signalp {signalp} && \\
Trinotate --db {trinotate_sqlite} --LOAD_infernal {infernal} && \\
Trinotate --db {trinotate_sqlite} --report -E {evalue} {other_options} \\
  > {trinotate_report}""".format(
            trinity_fasta=trinity_fasta,
            trinotate_sqlite=trinotate_sqlite,
            transdecoder_pep=transdecoder_pep,
            swissprot_blastx=swissprot_blastx,
            swissprot_blastp=swissprot_blastp,
            transdecoder_pfam=transdecoder_pfam,
            tmhmm=tmhmm,
            signalp=signalp,
            infernal=infernal,
            evalue=global_conf.global_get('trinotate', 'evalue'),
            other_options=global_conf.global_get('trinotate', 'other_options', required=False),
            trinotate_report=trinotate_report
        )
    )
