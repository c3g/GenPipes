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

def ensemble(input_callers, output, options):
    
    return Job(
        input_callers,
        [output],
        [
            ['bcbio_ensemble', 'module_bcbio_variation_recall'],
            ['bcbio_ensemble', 'module_bcftools'],
            ['bcbio_ensemble', 'module_java'],
        ],
        command="""\
$BCBIO_VARIATION_RECALL_HOME/bcbio.variation.recall ensemble \\
  {options} \\
  {output} \\
  {reference_sequence} \\
  {input_callers}""".format(
        tmp_dir=config.param('bcbio_ensemble', 'tmp_dir'),
        java_other_options=config.param('bcbio_ensemble', 'java_other_options'),
        ram=config.param('bcbio_ensemble', 'ram'),
        options=options,
        output=output if output else "-",
        reference_sequence=config.param('bcbio_ensemble', 'genome_fasta', param_type='filepath'),
        input_callers="  ".join("  \\\n  " + caller for caller in input_callers)
        )
    )

