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

def blastdbcmd(
    entry_file,
    entry_cmd,
    outfile
    ):

    return Job(
        [entry_file],
        [outfile],
        [
            ['blast_blastdbcmd', 'module_blast']
        ],
        command="""\
blastdbcmd \\
  -db {blast_db} \\
  -entry {entry_cmd} \\
  -outfmt %f \\
  > {outfile}""".format(
        blast_db=config.param('blast_blastdbcmd', 'blast_db', param_type='prefixpath'),
        entry_cmd=entry_cmd,
        outfile=outfile
    ))

def blastn_on_db(db, query, output, other_options=""):
    return Job(
        [query],
        [output],
        [['blast_blastnondb', 'module_blast']],
        command="""\
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

    return Job(
        [infile_fasta, coverage_bed],
        [outfile, os.path.join(outdir, "blastCov.tsv"), os.path.join(outdir, "contigsCoverage.tsv")],
        [
            ['blast_dcmegablast', 'module_blast'],
            ['blast_dcmegablast', 'module_R'],
            ['blast_dcmegablast', 'module_mugqic_tools']
        ],
        command="""\
blastn -task dc-megablast \\
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
        max_target_seqs=config.param('blast_dcmegablast', 'max_target_seqs', param_type='posint'),
        threads=config.param('blast_dcmegablast', 'threads', param_type='posint'),
        blast_db=config.param('blast_dcmegablast', 'blast_db', param_type='prefixpath'),
        outfile=outfile,
        coverage_bed=coverage_bed,
        outdir=outdir
    ))


# Parallel blast using fasta chunks
def parallel_blast(fasta, query, blast, program, db, cpu):    
    return(Job(
        [fasta],
        [blast],
        [['blast', 'module_perl'], ['blast', 'module_mugqic_tools'], ['blast', 'module_blast']],
        command="""\
parallelBlast.pl \\
-file {query} \\
--OUT {blast} \\
-n {cpu} \\
--BLAST "{program} -db {db} -max_target_seqs 1 -outfmt '6 std stitle'" """.format(
            query=query,
            blast=blast,
            cpu=cpu,
            program=program,
            db=db
        ))
    )
