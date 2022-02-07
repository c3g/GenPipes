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

def circularize(fasta_consensus, corrected_fastq, circlator_output):

    inputs = [fasta_consensus, corrected_fastq]
    outputs = [circlator_output]
    return Job(
        inputs,
        outputs,
        [
            ['circlator', 'module_python'],
            ['circlator', 'module_bwa'],
            ['circlator', 'module_samtools'],
            ['circlator', 'module_mummer'],
            ['circlator', 'module_spades'],
            ['circlator', 'module_prodigal'],
        ],
        command="""\
rm -rf {output} && \\
circlator all \\
  {fasta} \\
  {fastq} \\
  {output}""".format(
            fasta=fasta_consensus,
            fastq=corrected_fastq,
            output=circlator_output
        )
    )
