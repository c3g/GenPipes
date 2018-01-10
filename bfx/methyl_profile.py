#!/usr/bin/env python

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

# Python Standard Modules
import os

# MUGQIC Modules
from core.pipeline import *
from core.config import *
from core.job import *

def combine(input, output):
    return Job(
        [input],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_perl']
        ],
        command="""\
methylProfile.bismark.pl \\
 -i {input} \\
 -o {output}""".format(
            input=input,
            output=output
         )
    )

def cpg_stats(input, cg_stats, lambda_stats, puc19_stats):
    return Job(
        [input],
        [cg_stats, lambda_stats, puc19_stats],
        [
            ['DEFAULT', 'module_mugqic_tools']
        ],
        command="""\
bash cpgStats.sh \\
  {input} \\
  {cg_stats} \\
  {lambda_stats} \\
  {puc19_stats}""".format(
            input=input,
            cg_stats=cg_stats,
            lambda_stats=lambda_stats,
            puc19_stats=puc19_stats
        )
    )

def cpg_cov_stats(input, output):
    return Job(
        [input],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/CpG_coverageStats.py \\
 -i {input} \\
 -o {output}""".format(
            input=input,
            output=output
         )
    )

def metrics_report(sample_list, inputs, output, target_bed):
    if not isinstance(inputs, list):
        inputs=[inputs]
    
    return Job(
        inputs,
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools']
        ],
        command="""\
bash methylseq_metrics.sh \\
  {sample_list} \\
  {output_file} \\
  {targeted_flag}""".format(
            sample_list=",".join(sample_list),
            output_file=output,
            targeted_flag=1 if target_bed else 0
        )
    )

def ihec_metrics_report(sample_list, inputs, output, target_bed):
  if not isinstance(inputs, list):
        inputs=[inputs]
  
    return Job(
        inputs,
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools']
        ],
        command="""\
bash methylseq_metrics_for_ihec.sh \\
  {sample_list} \\
  {output_file} \\
  {targeted_flag}""".format(
            sample_list=",".join(sample_list),
            output_file=output,
            targeted_flag=1 if target_bed else 0
        )
    )
