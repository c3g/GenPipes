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
import logging

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)


def freebayes(input_bam, output_file, options=None, ini_section='freebayes'):
    # inputs = [input]
    # output = [prefix + ".vcf"]

    return Job(
        input_files=[input_bam],
        output_files=[output_file],
        module_entries=[
            [ini_section, 'module_freebayes']
        ],

        command="""\
freebayes -f {reference_genome} \\
  {bed_targets} \\
  {options} \\
  {other_options} \\
  {input_bam} \\
  | sed s/QR,Number=1,Type=Integer/QR,Number=1,Type=Float/ > {output_file}""".format(
    reference_genome=global_config_parser.param(ini_section, 'genome_fasta', param_type='filepath'),
    bed_targets="-t " + global_config_parser.param(ini_section, 'bed_targets', param_type='filepath', required=False) if global_config_parser.param(ini_section, 'bed_targets', required=False) else "",
    options=options if options else "",
    other_options=global_config_parser.param(ini_section, 'freebayes_options'),
    input_bam=input_bam,
    output_file=output_file
    )
        )


def process_gvcf(intput_gvcf, output_masks, output_ambiguous, output_consensus, ini_section='freebayes'):
    return Job(
        input_files=[intput_gvcf],
        output_files=[output_masks, output_ambiguous, output_consensus],
        module_entries=[
            [ini_section, 'module_freebayes'],
            [ini_section, 'module_python']
        ],

        command="""\
process_gvcf.py {other_options} -m {output_masks} -v {output_ambiguous} -c {output_consensus} {intput_gvcf}""".format(
            other_options=global_config_parser.param(ini_section, 'process_gvcf_options'),
            output_masks=output_masks,
            output_ambiguous=output_ambiguous,
            output_consensus=output_consensus,
            intput_gvcf=intput_gvcf
            )
        )
