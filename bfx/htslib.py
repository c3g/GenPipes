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

#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def bgzip_tabix_vcf(input, output):

    return Job(
        [input],
        [output],
        [
            ['DEFAULT', 'module_htslib'],
        ],
        command="""\
bgzip -cf \\
{input} > \\
{output} && tabix -pvcf {output} \\
        """.format(
        input=" \\\n " + input if input else "",
        output=output
        )
    )

def bgzip_tabix_bed(input, output):

    return Job(
        [input],
        [output],
        [
            ['DEFAULT', 'module_htslib'],
        ],
        command="""\
bgzip -cf \\
{input} \\
{output} && tabix -pbed {output} \\
        """.format(
        input=" \\\n  > " + input if input else "",
        output=output
        )
    )
