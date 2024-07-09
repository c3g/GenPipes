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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def run(tumor, patient_name, normal, bed, ini_section='svaba_run'):
    outputs = [patient_name + ".svaba.sv.vcf"]
    
    if normal:
        somatic_sv = patient_name + ".svaba.somatic.sv.vcf"
        germline_sv = patient_name + ".svaba.germline.sv.vcf"
        outputs = [somatic_sv, germline_sv]

    return Job(
        [tumor, normal],
        outputs,
        [
            [ini_section, 'module_svaba']
        ],
        command="""\
svaba run {options} \\
        -G {ref} \\
        -D {dbsnp}{bed} \\
        -a {name} \\
        -t {tumor} \\
        {normal}""".format(
            options=global_conf.global_get(ini_section, 'options'),
            ref=global_conf.global_get(ini_section, 'ref', param_type='filepath'),
            dbsnp=global_conf.global_get(ini_section, 'dbsnp'),
            bed=" -k " + bed if bed else "",
            name=patient_name,
            normal="-n " + normal if normal else "",
            tumor=tumor,
        )
    )
