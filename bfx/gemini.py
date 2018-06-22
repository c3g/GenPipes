################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
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

#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
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
        options=config.param('gemini_annotations', 'options'),
        variants=variants,
        output=gemini_output,
        temp=tmp_dir
        )
    )

def vcfanno(input, output):

    return Job(
        [input],
        [output],
        [
            ['gemini_annotations', 'module_vcfanno'],
            ['gemini_annotations', 'module_tabix']
        ],
        command="""\
vcfanno {options} \\
-base-path {gemini_data} \\
-lua {rare_disease_lua} {rare_disease_GRCh38_conf} \\
{input} | \\
bgzip -c > {output}""".format(
        options=config.param('gemini_annotations', 'vcfanno_options'),
        gemini_data=config.param('gemini_annotations', 'gemini_data'),
        rare_disease_lua=config.param('gemini_annotations', 'rare_disease_lua'),
        rare_disease_GRCh38_conf=config.param('gemini_annotations', 'rare_disease_GRCh38_conf'),
        input=input,
        output=output
        )
    )

def vcf2db(input, output_db):

    return Job(
        [input],
        [output_db],
        [['gemini_annotations', 'module_python']],
        command="""\
python {vcf2db_script} --legacy-compression \\
{input} {ped_file} \\
{output}""".format(
        vcf2db_script=config.param('gemini_annotations', 'vcf2db_script'),
        input=input,
        ped_file=config.param('gemini_annotations', 'ped_file'),
        output=output_db
        )
    )
