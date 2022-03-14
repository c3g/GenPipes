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

def scalpel_somatic(inputNormal, inputTumor, outputDir, bed):

    return Job(
        [inputNormal, inputTumor, bed],
        [os.path.join(outputDir, 'main', 'somatic.5x.indel.vcf'), os.path.join(outputDir, 'main', 'common.5x.indel.vcf')],
        [
            ['scalpel', 'module_perl'],
            ['scalpel', 'module_scalpel']
        ],
        command="""\
scalpel-discovery --somatic \\
  --ref {reference_sequence} \\
  --normal {inputNormal} \\
  --tumor {inputTumor} \\
  --dir {outputDir} \\
  --numprocs {cores_per_job} \\
  --bed {bed}""".format(
        reference_sequence=global_config_parser.param('scalpel', 'genome_fasta', param_type='filepath'),
        inputNormal=inputNormal,
        inputTumor=inputTumor,
        outputDir=outputDir,
        cores_per_job=global_config_parser.param('scalpel', 'cores_per_job'),
        bed=bed
        )
    )

def scalpel_somatic_2pass(inputNormal, inputTumor, outputDir, bed):

    return Job(
        [inputNormal, inputTumor, bed],
        [os.path.join(outputDir, 'main', 'somatic.5x.indel.vcf'), os.path.join(outputDir, 'main', 'common.5x.indel.vcf')],
        [
            ['scalpel', 'module_perl'],
            ['scalpel', 'module_scalpel']
        ],
        command="""\
scalpel-discovery --somatic --two-pass \\
  --ref {reference_sequence} \\
  --normal {inputNormal} \\
  --tumor {inputTumor} \\
  --dir {outputDir} \\
  --numprocs {cores_per_job} \\
  --bed {bed}""".format(
        reference_sequence=global_config_parser.param('scalpel', 'genome_fasta', param_type='filepath'),
        inputNormal=inputNormal,
        inputTumor=inputTumor,
        outputDir=outputDir,
        cores_per_job=global_config_parser.param('scalpel', 'cores_per_job'),
        bed=bed
        )
    )
