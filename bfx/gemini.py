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
from core.job import *


def gemini_annotations(variants, gemini_output, tmp_dir):
    return Job(
        [variants],
        [gemini_output],
        [
            ['gemini_annotations', 'module_gemini'],
            ['gemini_annotations', 'module_htslib']
        ],
        command="""\
gemini load -v {variants} \\
  {options} \\
  --tempdir {temp} \\
  {output}""".format(
            options=global_config_parser.param('gemini_annotations', 'options'),
            variants=variants,
            output=gemini_output,
            temp=tmp_dir
        )
    )


def set_somatic(ped, database, output):
    return Job(
        [database],
        [output],
        [
            ['gemini_annotations', 'module_gemini'],
            ['gemini_annotations', 'module_htslib']
        ],
        command="""\
gemini amend \\
  --sample {ped} \\
  {database} && \\
gemini set_somatic \\
  {options} \\
  {database} > \\
  {output}""".format(
            options=global_config_parser.param('set_somatic_and_actionable_mutations', 'set_somatic'),
            ped=ped,
            database=database,
            output=output,
        )
    )


def actionable_mutations(database, output):
    return Job(
        [database],
        [output],
        [
            ['gemini_annotations', 'module_gemini'],
            ['gemini_annotations', 'module_htslib']
        ],
        command="""\
gemini actionable_mutations \\
  {database} \\
  > {output}""".format(
            database=database,
            output=output,
        )
    )
