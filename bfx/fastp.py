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
import configparser
import logging
import os
import re

# MUGQIC Modules
from core.config import config
from core.job import Job, concat_jobs

def basic_qc(
    input1,
    input2,
    output_json_path,
    output_html_path=None,
    overrepresentation_analysis=True
    ):

    if input2:
        inputs = [input1, input2]
    else:
        inputs = [input1]

    num_threads = config.param('fastp', 'threads', required=False, param_type='posint')
    output_files = filter(None, [output_json_path, output_html_path])

    return Job(
        inputs,
        output_files,
        module_entries=[
            ['fastp', 'module_fastp']
        ],
        command="""\
fastp -V \\
  {rds1} \\
  {rds2} \\
  {cpus} \\
  {json} \\
  {html}""".format(
            rds1 = "\\\n  --in1 "    + input1,
            rds2 = "\\\n  --in2 "    + input2 if input2 else "",
            orep = "\\\n  --overrepresentation_analysis" if overrepresentation_analysis else "",
            cpus = "\\\n  --thread " + str(num_threads) if num_threads else "",
            json = "\\\n  --json "   + output_json_path if output_json_path else "",
            html = "\\\n  --html "   + output_html_path if output_html_path else "",
        )
    )

def trim(
        adapter,
        input1,
        input2,
        output1,
        output2,
        output_json_path,
        output_html_path=None,
        ini_section = 'DEFAULT'
):

    if input2:
        inputs = [input1, input2]
    else:
        inputs = [input1]
        
    if output2:
        outputs = [output1, output2, output_json_path, output_html_path]
    else:
        outputs = [output1, output_json_path, output_html_path]

    num_threads = config.param(ini_section, 'threads', required=False, param_type='posint')
    paramters = config.param(ini_section, 'options', required=False)

    return Job(
        inputs,
        outputs,
        module_entries=[
            ['fastp', 'module_fastp']
        ],
        command="""\
fastp -V \\
  {options} \\
  {adapter} \\
  {rds1} \\
  {rds2} \\
  {out1} \\
  {out2} \\
  {cpus} \\
  {json} \\
  {html}""".format(
            options = "\\\n  "                  + paramters,
            adapter = "\\\n  --adapter_fasta "  + adapter,
            rds1 = "\\\n  --in1 "               + input1,
            rds2 = "\\\n  --in2 "               + input2 if input2 else "",
            out1 = "\\\n  --out1 "              + output1,
            out2 = "\\\n  --out2 "              + output2 if output2 else "",
            cpus = "\\\n  --thread "            + str(num_threads) if num_threads else "",
            json = "\\\n  --json "              + output_json_path if output_json_path else "",
            html = "\\\n  --html "              + output_html_path if output_html_path else "",
        )
    )