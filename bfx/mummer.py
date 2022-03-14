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

    job = Job(
        [fasta_reference, fasta_consensus],
        [outfile, outfile + ".png", outfile2],
        [
            ['mummer_reference', 'module_perl'],
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
        c=global_config_parser.param('mummer_reference', 'c', param_type='posint'),
        prefix1=prefix1,
        fasta_reference=fasta_reference,
        fasta_consensus=fasta_consensus,
        title=title,
        prefix2=prefix2,
        outfile=outfile,
        prefix3=prefix3,
        x=global_config_parser.param('mummer_reference', 'x', param_type='posint'),
        infile2=infile2,
        outfile2=outfile2
    ))

    # Mammouth does not have libgd by default. Module must be loaded explicitely
    if global_config_parser.param('mummer_reference', 'module_libgd', required=False):
        job.modules.append(global_config_parser.param('mummer_reference', 'module_libgd', required=False))

    return job

def self(
    prefix1,
    fasta_consensus,
    title,
    prefix2,
    outfile
    ):

    job = Job(
        [fasta_consensus],
        [outfile, outfile + ".png"],
        [
            ['mummer_self', 'module_perl'],
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
        c=global_config_parser.param('mummer_self', 'c', param_type='posint'),
        prefix1=prefix1,
        fasta_consensus=fasta_consensus,
        title=title,
        prefix2=prefix2,
        outfile=outfile
    ))

    # Mammouth does not have libgd by default. Module must be loaded explicitely
    if global_config_parser.param('mummer_self', 'module_libgd', required=False):
        job.modules.append(global_config_parser.param('mummer_self', 'module_libgd', required=False))

    return job
