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

# MUGQIC Modules
from core.config import config
from core.job import Job

def fastp_basic_qc(input1, input2, output_json_path, output_html_path=None, overrepresentation_analysis=True):
    if input2:
        inputs = [input1, input2]
    else:
        inputs = [input1]

    output_files = [output_json_path, output_html_path]

    return Job(
        inputs,
        output_files,
        module_entries=[['fastp', 'module_fastp']],
        command="""\
fastp --in1 {rds1}{rds2}{orep} \\
    --thread {cpus} \\
    --json {json}{html}""".format(
    rds1=input1,
    rds2=" --in2 " + input2 if input2 else "",
    orep=" --overrepresentation_analysis " if overrepresentation_analysis else "",
    cpus=config.param('fastp', 'threads', required=False, type='posint'),
    json=output_json_path,
    html=" --html " + output_html_path if output_html_path else ""
    )
        )
