################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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


def report(
    input,
    output_dir,
    tumor_id,
    ini_section ='report_cpsr'
    ):

    assembly = global_conf.global_get(ini_section, 'assembly')
    output = [
        os.path.join(output_dir, tumor_id + ".cpsr." + assembly + ".json.gz"),
        os.path.join(output_dir, tumor_id + ".cpsr." + assembly + ".html")
    ]

    if global_conf.global_get(ini_section, 'module_pcgr').split("/")[2] >= "1":
        call = 'cpsr'
        module = 'module_pcgr'
    else:
        call = 'cpsr.py'
        module = 'module_cpsr'
    if global_conf.global_get(ini_section, 'module_pcgr').split("/")[2] >= "2":
        return report2(input, output_dir, tumor_id, ini_section)
    else:
        return Job(
            [input],
            output,
            [
                [ini_section, module],
            ],
            command="""\
{call} {options} \\
    --input_vcf {input} \\
    --pcgr_dir $PCGR_DATA \\
    --output_dir {output_dir} \\
    --genome_assembly {assembly} \\
    --sample_id {tumor_id}""".format(
            call=call,
            options=global_conf.global_get(ini_section, 'options'),
            input=input,
            output_dir=output_dir,
            assembly=global_conf.global_get(ini_section, 'assembly'),
            tumor_id=tumor_id
        )
    )

def report2(
        input, 
        output_dir, 
        tumor_id, 
        ini_section
        ):
    assembly = global_conf.global_get(ini_section, 'assembly')
    output = [
        os.path.join(output_dir, tumor_id + ".cpsr." + assembly + ".vcf.gz"),
        os.path.join(output_dir, tumor_id + ".cpsr." + assembly + ".html"),
        os.path.join(output_dir, tumor_id + ".cpsr." + assembly + ".classification.tsv.gz"),
        os.path.join(output_dir, tumor_id + ".cpsr." + assembly + ".conf.yaml")
    ]
    # use tmp dir for cpsr to avoid disk quota issues caused by bcftools tmp dir settings
    return Job(
        [input],
        output,
        [
            [ini_section, 'module_pcgr'],
        ],
        command="""\
mkdir -p {tmp_dir}/cpsr && \\
cpsr {options} \\
    --input_vcf {input} \\
    --refdata_dir $PCGR_DATA \\
    --vep_dir $PCGR_VEP_CACHE \\
    --output_dir {tmp_dir}/cpsr \\
    --genome_assembly {assembly} \\
    --sample_id {tumor_id} && \\
cp -r {tmp_dir}/cpsr {output_dir}""".format(
            options=global_conf.global_get(ini_section, 'options_v2'),
            input=input,
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            output_dir=os.path.dirname(output_dir),
            assembly=global_conf.global_get(ini_section, 'assembly'),
            tumor_id=tumor_id
        )
    )
