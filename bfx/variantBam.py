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

def run(input, output, region=None, exclude_region=None):
    return Job(
        [input],
        [output],
        [
            ['variantBam', 'module_variantBam'],
        ],
        command="""\
variant \\
  {input} \\
  {options}{region}{exclude_region} \\
  --reference {reference_fasta} \\
  --output {output}""".format(
        input=input,
        options=global_config_parser.param('samtools_cram_output', 'variantBam_options'),
        reference_fasta=global_config_parser.param('samtools_cram_output', 'genome_fasta', param_type='filepath'),
        region=" \\\n  " + region if region else "",
        exclude_region=" \\\n  " + exclude_region if exclude_region else "",
        output=output
        )
    )
