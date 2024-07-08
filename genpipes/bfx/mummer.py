################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def reference(
    prefix1,
    fasta_reference,
    fasta_consensus,
    title,
    prefix2,
    outfile,
    prefix3,
    infile2,
    outfile2,
    ini_section='mummer_reference'
    ):

    job = Job(
        [fasta_reference, fasta_consensus],
        [outfile, outfile + ".png", outfile2],
        [
            [ini_section, 'module_perl'],
            [ini_section, 'module_mummer'],
            [ini_section, 'module_gnuplot']
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
        c=global_conf.global_get(ini_section, 'c', param_type='posint'),
        prefix1=prefix1,
        fasta_reference=fasta_reference,
        fasta_consensus=fasta_consensus,
        title=title,
        prefix2=prefix2,
        outfile=outfile,
        prefix3=prefix3,
        x=global_conf.global_get(ini_section, 'x', param_type='posint'),
        infile2=infile2,
        outfile2=outfile2
    ))

    # Mammouth does not have libgd by default. Module must be loaded explicitely
    if global_conf.global_get(ini_section, 'module_libgd', required=False):
        job.modules.append(global_conf.global_get(ini_section, 'module_libgd', required=False))

    return job

def self(
    prefix1,
    fasta_consensus,
    title,
    prefix2,
    outfile,
    ini_section='mummer_self'
    ):

    job = Job(
        [fasta_consensus],
        [outfile, outfile + ".png"],
        [
            [ini_section, 'module_perl'],
            [ini_section, 'module_mummer'],
            [ini_section, 'module_gnuplot']
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
        c=global_conf.global_get(ini_section, 'c', param_type='posint'),
        prefix1=prefix1,
        fasta_consensus=fasta_consensus,
        title=title,
        prefix2=prefix2,
        outfile=outfile
    ))

    # Mammouth does not have libgd by default. Module must be loaded explicitely
    if global_conf.global_get(ini_section, 'module_libgd', required=False):
        job.modules.append(global_conf.global_get(ini_section, 'module_libgd', required=False))

    return job
