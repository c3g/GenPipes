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
        tmp_dir=global_conf.global_get('bcbio_ensemble', 'tmp_dir'),
        java_other_options=global_conf.global_get('bcbio_ensemble', 'java_other_options'),
        ram=global_conf.global_get('bcbio_ensemble', 'ram'),
        options=options,
        output=output if output else "-",
        reference_sequence=global_conf.global_get('bcbio_ensemble', 'genome_fasta', param_type='filepath'),
        input_callers="  ".join("  \\\n  " + caller for caller in input_callers)
        )
    )

