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
from core.config import *
from core.job import *

def run(tumor, patient_name, normal, bed):
    outputs = [patient_name + ".svaba.sv.vcf"]
    
    if normal:
        somatic_sv = patient_name + ".svaba.somatic.sv.vcf"
        germline_sv = patient_name + ".svaba.germline.sv.vcf"
        outputs = [somatic_sv, germline_sv]

    return Job(
        [tumor, normal],
        outputs,
        [
            ['svaba_run', 'module_svaba'],
#            ['svaba_run', 'module_gcc']
        ],
        command="""\
svaba run {options} \\
        -G {ref} \\
        -D {dbsnp}{bed} \\
        -a {name} \\
        -t {tumor} \\
        {normal}""".format(
            options=global_config_parser.param('svaba_run', 'options'),
            ref=global_config_parser.param('svaba_run', 'ref', param_type='filepath'),
            dbsnp=global_config_parser.param('svaba_run', 'dbsnp'),
            bed=" -k " + bed if bed else "",
            name=patient_name,
            normal="-n " + normal if normal else "",
            tumor=tumor,
        )
    )
