#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def blastdbcmd(
    entry_file,
    entry_cmd,
    outfile
    ):

    return  Job(
        [entry_file],
        [outfile],
        [
            ['blast_blastdbcmd', 'module_memtime'],
            ['blast_blastdbcmd', 'module_blast']
        ],
        command = """\
memtime blastdbcmd \\
  -db {blast_db} \\
  -entry {entry_cmd} \\
  -outfmt %f \\
  > {outfile}""".format(
        blast_db=config.param('blast_blastdbcmd', 'blast_db', type='prefixpath'),
        entry_cmd=entry_cmd,
        outfile=outfile
    ))

def blastn_on_db(db, query, output, other_options=""):
    return Job(
        [query],
        [output],
        [['DEFAULT', 'module_blast']],
        command = """\
blastn {other_options} \\
  -db {db} \\
  -query {query} \\
  -out {output}""".format(
        other_options=other_options,
        db=db,
        query=query,
        output=output
    ))

def dcmegablast(
    infile_fasta,
    outfmt,
    outfile,
    coverage_bed,
    outdir
    ):

    tmp_outfile = os.path.splitext(outfile)[0] + ".all.tmp"

    return  Job(
        [infile_fasta],
        [outfile],
        [
            ['blast_dcmegablast', 'module_memtime'],
            ['blast_dcmegablast', 'module_blast'],
            ['blast_dcmegablast', 'module_R'],
            ['blast_dcmegablast', 'module_mugqic_tools']
        ],
        command = """\
memtime blastn -task dc-megablast \\
  -query {infile_fasta} \\
  -outfmt "{outfmt} qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle sskingdoms sscinames scomnames" \\
  -out {tmp_outfile} \\
  -max_target_seqs {max_target_seqs} \\
  -num_threads {threads} \\
  -db {blast_db} && \\
pacBioKeepBlastBestHits.pl \\
  --infile {tmp_outfile} \\
  --n {max_target_seqs} \\
  > {outfile} && \\
pacBioMergeCovToBlast.R \\
  -c {coverage_bed} \\
  -b {outfile} \\
  -o {outdir}""".format(
        infile_fasta=infile_fasta,
        outfmt=outfmt,
        tmp_outfile=tmp_outfile,
        max_target_seqs=config.param('blast_dcmegablast', 'max_target_seqs', type='posint'),
        threads=config.param('blast_dcmegablast', 'threads', type='posint'),
        blast_db=config.param('blast_dcmegablast', 'blast_db', type='prefixpath'),
        outfile=outfile,
        coverage_bed=coverage_bed,
        outdir=outdir
    ))
