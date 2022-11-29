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
import logging

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

### Start from here ###

def bamcoverage(input_bam, output_file):
    return Job(
        input_files=[input_bam],
        output_files=[output_file],
        module_entries=[['wiggle', 'module_deeptools'] # need to add deeptools to the .ini***
        ],
        # name=job_name,
        command="""\
bamCoverage --verbose \\
    --outFileFormat bigwig \\
    {other_options} \\
    --bam ${input_bam} \\
    --outFileName ${output_file}.bw """.format(
        output_file=output_file,
        input_bam=input_bam,
        other_options=config.param('wiggle', 'other_options', required=True),
        )
  )
