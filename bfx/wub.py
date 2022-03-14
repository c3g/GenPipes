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


def bam_alignment_qc(bam_file, pickle, ini_section="wub_metrics"):
    """
    """

    return Job(
        [bam_file],
        [pickle],
        [[ini_section, "module_wub"]],
        command="""\
bam_alignment_qc.py -f {reference} \\
  -p {pickle} \\
  {bam_file}""".format(
            reference=global_config_parser.param(ini_section, 'genome_fasta', required=True),
            pickle=pickle,
            bam_file=bam_file
        )
    )
