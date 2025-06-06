################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
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
import logging

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

log = logging.getLogger(__name__)


def trim_primers(input_bam, prefix):
    inputs = [input_bam]
    output = [prefix + ".bam"]

    return Job(
        input_files=inputs,
        output_files=output,
        module_entries=[
            ['ivar', 'module_ivar'],
            ['htslib', 'module_htslib']
        ],

        command="""\
ivar trim -i {input_bam} \\
  -p {prefix} \\
  {bed_file} \\
  {primer_pair} \\
  {other_options}""".format(
            input_bam=input_bam,
            prefix=prefix,
            bed_file="-b " + global_conf.global_get('ivar_trim_primers', 'bed_primers', param_type='filepath') if global_conf.global_get('ivar_trim_primers', 'bed_primers') else "",
            primer_pair="-f " + global_conf.global_get('ivar_trim_primers', 'tsv_primer_pair', param_type='filepath') if global_conf.global_get('ivar_trim_primers', 'tsv_primer_pair') else "",
            other_options=global_conf.global_get('ivar_trim_primers', 'other_options')
            ),
        )

def call_variants(prefix):
    # inputs = [input]
    output = [prefix + ".tsv"]

    return Job(
        input_files=[],
        output_files=output,
        module_entries=[
            ['ivar', 'module_ivar'],
            ['htslib', 'module_htslib']
        ],

        command="""\
ivar variants -p {prefix} \\
  -r {reference_genome} \\
  {gff_file} \\
  {other_options}""".format(
            prefix=prefix,
            reference_genome=global_conf.global_get('ivar_call_variants', 'genome_fasta', param_type='filepath'),
            gff_file="-g " + global_conf.global_get('ivar_call_variants', 'gff_orf', param_type='filepath') if global_conf.global_get('ivar_call_variants', 'gff_orf') else "",
            other_options=global_conf.global_get('ivar_call_variants', 'other_options')
            ),
        )

def create_consensus(prefix):
    # inputs = [input]
    output = [prefix + ".fa", prefix + ".qual.txt"]

    return Job(
        input_files=[],
        output_files=output,
        module_entries=[
            ['ivar', 'module_ivar'],
            ['htslib', 'module_htslib']
        ],

        command="""\
ivar consensus -p {prefix} {other_options}""".format(
            prefix=prefix,
            other_options=global_conf.global_get('ivar_create_consensus', 'other_options')
            ),
        )

def tsv_to_vcf(input, output):
    input = [input]
    output = [output]

    return Job(
        input_files=input,
        output_files=output,
        module_entries=[
            ['ivar', 'module_ivar'],
            ['htslib', 'module_htslib']
        ],

        command="""\
ivar_variants_to_vcf.py {input} {output}""".format(
            input=" \\\n  ".join(input),
            output=" \\\n  ".join(output)
            ),
        )
