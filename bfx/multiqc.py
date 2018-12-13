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

# MUGQIC Modules
from core.config import *
from core.job import *


def mutliqc_run(yamlFile, input_files):


    command = """export MULTIQC_CONFIG_PATH={yamlFile} && \\
    multiqc .""".format(
        yamlFile = yamlFile)
    ## for now multiqc will run after hicup alignments are complete. Once Homer is added to mutliqc, the input must change to refect homer tag dirs

    return Job(input_files = input_files,
            output_files = ["Analysis_Summary_Report.html"],
            module_entries = [["multiqc_report", "module_multiqc"]],
            name = "multiqc_report",
            command = command
            )

def run(inputs, output, input_dep=False):
    return Job(
        input_dep,
        [output],
        [
            ['multiqc', 'module_multiqc'],
        ],
        command="""\
multiqc -f {options} \\
{input} \\
-n {output}""".format(
            options=config.param('multiqc', 'options'),
            input=" ".join([" \\\n  " + input for input in inputs]),
            output=output,
            )
        )