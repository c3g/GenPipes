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

def somatic(input, normal_name, tumor_name, output):
    return Job(
        [input],
        [output],
        [
            ['vawk', 'module_vawk'],
        ],
        command="""\
zcat {input} | \\ 
        vawk --header \\
        '(S${tumor_name}$GT!="0/0" && S${tumor_name}$GT!="./." \\
        && S${tumor_name}$GT!=S${normal_name}$GT) \\
        && (S${normal_name}$GT=="0/0" || S${normal_name}$GT=="./.")' \\
        {output}""".format(
            input=input,
            normal_name=normal_name,
            tumor_name=tumor_name,
            output="> " + output if output else "",
        )
    )

def germline(input, normal_name, tumor_name, output):
    return Job(
        [input],
        [output],
        [
            ['vawk', 'module_vawk'],
        ],
        command="""\
zcat {input} | \\ 
        vawk --header \\
        '(S${normal_name}$GT!="0/0" && S${normal_name}$GT!="./." && S${tumor_name}$GT!="./.")' \\
        {output}""".format(
            input=input,
            normal_name=normal_name,
            tumor_name=tumor_name,
            output="> " + output if output else "",
        )
    )