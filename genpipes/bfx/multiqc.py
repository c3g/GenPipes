################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
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
import os

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def run(inputs, output, ini_section='multiqc'):
    if not isinstance(inputs, list):
        inputs = [inputs]
    output = output + ".html"
    return Job(
        inputs,
        [output],
        [
            [ini_section, 'module_multiqc']
        ],
        command="""\
mkdir -p {tmp_output} && \\
TMPDIR={tmp_output} && \\
multiqc -f {options} \\
{input} \\
-n {output} && \\
rm -rf {tmp_output}""".format(
            options=global_conf.global_get(ini_section, 'options', required=False),
            input=" ".join([" \\\n  " + input for input in inputs]),
            output=output,
            tmp_output=os.path.join(os.path.dirname(output), "tmp_multiqc")
        )
    )

def multiqc_run(
    yamlFile,
    input_files
    ):
    ## for now multiqc will run after hicup alignments are complete. Once Homer is added to multiqc, the input must change to refect homer tag dirs
    return Job(
        input_files,
        ["Analysis_Summary_Report.html"],
        module_entries = [
            ["multiqc_report", "module_multiqc"]
        ],
        command = """\
export MULTIQC_CONFIG_PATH={yamlFile} && \\
  multiqc .""".format(
            yamlFile = yamlFile
        )
    )
