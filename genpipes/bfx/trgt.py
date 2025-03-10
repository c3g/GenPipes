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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job


def genotype(
    input_file,
    repeat_file,
    output_prefix,
    ini_section='trgt_genotyping'
    ):
    """
    Call tandem repeats with TRGT.

    :return: a job for trgt genotyping
    """

    genome_fasta = global_conf.global_get(ini_section, 'genome_fasta', required=True)
    outputs = [
        output_prefix + ".vcf.gz",
        output_prefix + ".spanning.bam"
    ]

    return Job(
        [input_file],
        outputs,
        [[ini_section, "module_trgt"]],
        command="""\
trgt genotype {other_options} \\
  --genome {genome_fasta} \\
  --reads {input_file} \\
  --repeats {repeat_file} \\
  --threads {threads} \\
  --output-prefix {output_prefix}""".format(
            threads=global_conf.global_get(ini_section, 'threads'),
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            genome_fasta=genome_fasta,
            input_file=input_file,
            repeat_file=repeat_file,
            output_prefix=output_prefix
        ),
        removable_files=outputs
    )
