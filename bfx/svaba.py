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

def run(tumor, patient_name, normal=None):
    somatic_sv = patient_name + ".svaba.somatic.sv.vcf"
    germline_sv = patient_name + ".svaba.germline.sv.vcf"
    output = [somatic_sv, germline_sv]

    return Job(
        [normal, tumor],
        output,
        [
            ['svaba_run', 'module_svaba'],
            ['svaba_run', 'module_gcc']
        ],
        command="""\
svaba run {options} \\
        -G {ref} \\
        {dbsnp} \\
        -a {name} \\
        -t {tumor} \\
        {normal}""".format(
            options=config.param('svaba_run', 'options'),
            ref=config.param('svaba_run', 'ref', type='filepath'),
            dbsnp=config.param('svaba_run', 'dbsnp'),
            name=patient_name,
            normal="-n " + normal if normal else "",
            tumor=tumor,
        )
    )