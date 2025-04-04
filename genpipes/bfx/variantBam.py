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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

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
        options=global_conf.global_get('samtools_cram_output', 'variantBam_options'),
        reference_fasta=global_conf.global_get('samtools_cram_output', 'genome_fasta', param_type='filepath'),
        region=" \\\n  " + region if region else "",
        exclude_region=" \\\n  " + exclude_region if exclude_region else "",
        output=output
        )
    )
