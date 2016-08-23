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

#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def decompose_and_normalize_mnps(input, vt_output):

    return Job(
        [input],
        [vt_output],
        [
            ['DEFAULT', 'module_htslib'],
            ['decompose_and_normalize_mnps', 'module_vt']
        ],
        command="""\
zcat {input} | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r {reference_sequence} - | bgzip -cf > {vt_output} \\
        """.format(
        input=input,
        reference_sequence=config.param('DEFAULT', 'genome_fasta', type='filepath'),
        vt_output=vt_output
        )
    )
