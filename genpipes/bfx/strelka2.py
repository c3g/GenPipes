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
# along with GenPipes. If not, see <http://www.gnu.org/licenses/>.
################################################################################

import os
import re

from ..core.config import global_conf
from ..core.job import Job

def somatic_config(
    input_normal,
    input_tumor,
    output_dir,
    callRegions=None,
    mantaIndels=None,
    ini_section='strelka2_paired_somatic'
    ):

    inputs = [input_normal, input_tumor]
    if callRegions:
        inputs.append(callRegions)
    if mantaIndels:
        inputs.append(mantaIndels)
    return Job(
        inputs,
        [os.path.join(output_dir, "runWorkflow.py")],
        [
            [ini_section, 'module_python'],
            [ini_section, 'module_strelka2']
        ],
        command="""\
rm -f {workflow} && \\
python $STRELKA2_HOME/bin/configureStrelkaSomaticWorkflow.py \\
  --normalBam {normal} \\
  --tumorBam {tumor} \\
  --referenceFasta {genome} \\
  {experiment_type} {callRegions} {mantaIndels} \\
  --runDir {output}""".format(
            workflow=os.path.join(output_dir, "runWorkflow.py"),
            normal=input_normal,
            tumor=input_tumor,
            genome=global_conf.global_get(ini_section,'genome_fasta', param_type='filepath'),
            experiment_type=global_conf.global_get(ini_section,'experiment_type_option', required=False),
            callRegions="\\\n  --callRegions " + callRegions if callRegions else "",
            mantaIndels="\\\n  --indelCandidates " + mantaIndels if mantaIndels else "",
            output=output_dir
        )
    )

def germline_config(
    input_normal,
    output_dir,
    callRegions=None,
    ini_section='strelka2_paired_germline'
    ):

    if not isinstance(input_normal, list):
        input_normal = [input_normal]

    return Job(
        input_normal,
        [os.path.join(output_dir, "runWorkflow.py")],
        [
            [ini_section, 'module_python'],
            [ini_section, 'module_strelka2']
        ],
        command="""\
rm -f {workflow} && \\
python $STRELKA2_HOME/bin/configureStrelkaGermlineWorkflow.py \\
  {normal} \\
  --referenceFasta {genome} \\
  {experiment_type} {callRegions} \\
  --runDir {output}""".format(
            workflow=os.path.join(output_dir, "runWorkflow.py"),
            normal="".join(" \\\n  --bam " + input for input in input_normal),
            genome=global_conf.global_get(ini_section,'genome_fasta', param_type='filepath'),
            experiment_type=global_conf.global_get(ini_section,'experiment_type_option', required=False),
            callRegions="\\\n  --callRegions " + callRegions if callRegions else "",
            output=output_dir
        )
    )

def run(
    input_dir,
    output_dep=[],
    ini_section='strelka2_paired_somatic'
    ):

    ram = global_conf.global_get(ini_section, 'ram')
    ram_num = re.match('[0-9]+', ram)
    ram_GB = ram_num.group()
    if 'm' in ram.lower():
        ram_GB = ram_num / 1024
    elif 't' in ram.lower():
        ram_GB = ram_num * 1024


    return Job(
        [os.path.join(input_dir, "runWorkflow.py")],
        output_dep,
        [
            [ini_section, 'module_python'],
            [ini_section, 'module_strelka2']
        ],
        command="""\
python {input_dir}/runWorkflow.py \\
  -m {mode}  \\
  -j {nodes} \\
  -g {ram} \\
  --quiet""".format(
            input_dir=input_dir,
            mode=global_conf.global_get(ini_section,'option_mode'),
            nodes=global_conf.global_get(ini_section,'option_nodes'),
            ram=ram_GB
        )
    )
