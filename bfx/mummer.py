#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def reference(
    prefix1,
    fasta_reference,
    fasta_consensus,
    title,
    prefix2,
    outfile,
    prefix3,
    infile2,
    outfile2
    ):

    return Job(
        [fasta_reference, fasta_consensus],
        [outfile, outfile + ".png", outfile2],
        [
            ['mummer_reference', 'module_mummer'],
            ['mummer_reference', 'module_gnuplot']
        ],
        command="""\
promer --maxmatch \\
  --mincluster {c} \\
  --prefix {prefix1} \\
  {fasta_reference} \\
  {fasta_consensus} && \\
mummerplot --png --layout --filter \\
  --title {title} \\
  --prefix {prefix2} \\
  {outfile} && \\
dnadiff \\
  --prefix {prefix3} \\
  {fasta_reference} \\
  {fasta_consensus} && \\
show-snps -rlTC \\
  -x {x} \\
  {infile2} \\
  > {outfile2}""".format(
        c=config.param('mummer', 'c', type='posint'),
        prefix1=prefix1,
        fasta_reference=fasta_reference,
        fasta_consensus=fasta_consensus,
        title=title,
        prefix2=prefix2,
        outfile=outfile,
        prefix3=prefix3,
        x=config.param('mummer', 'x', type='posint'),
        infile2=infile2,
        outfile2=outfile2
    ))

def self(
    prefix1,
    fasta_consensus,
    title,
    prefix2,
    outfile
    ):

    return Job(
        [fasta_consensus],
        [outfile, outfile + ".png"],
        [
            ['mummer_self', 'module_mummer'],
            ['mummer_self', 'module_gnuplot']
        ],
        command="""\
nucmer --maxmatch \\
  --mincluster {c} \\
  --prefix {prefix1} \\
  {fasta_consensus} \\
  {fasta_consensus} && \\
mummerplot --png --layout --filter \\
  --title {title} \\
  --prefix {prefix2} \\
  {outfile}""".format(
        c=config.param('mummer', 'c', type='posint'),
        prefix1=prefix1,
        fasta_consensus=fasta_consensus,
        title=title,
        prefix2=prefix2,
        outfile=outfile
    ))
