################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
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


def run(input, output, ini_section='run_vcfanno'):
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_vcfanno'],
        ],

        command="""\\
vcfanno {options} \\
    -lua {lua} {config} \\
    {input} \\
    {output}""".format(
            options=global_conf.global_get(ini_section, 'options'),
            lua=global_conf.global_get(ini_section, 'lua'),
            config=global_conf.global_get(ini_section, 'config'),
            input=input,
            output=" \\\n > " + output if output else "",
        )
    )
