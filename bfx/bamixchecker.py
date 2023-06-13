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
def run_bam(
    inputfiles,
    outputfiles,
    out_dir,
    inputlist,
    reference,
    options,
    ini_section='run_bamixchecker'
    ):
    return Job(
        inputfiles,
        outputfiles,
        [
            [ini_section, 'module_python'],
            [ini_section, 'module_java'],
            [ini_section, 'module_bedtools'],
            [ini_section, 'module_R'],
            [ini_section, 'module_pandoc'],
            [ini_section, 'module_gatk'],
            [ini_section, 'module_bamixchecker']
        ],
        command="""\
python $BAMixChecker_PATH/BAMixChecker.py \\
  -l {inputlist} \\
  -r {reference} \\
  -o {out_dir} \\
  {options}""".format(
            inputlist=inputlist,
            out_dir=out_dir,
            options=options,
            reference=reference
        )
    )