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
import os

# MUGQIC Modules
from core.job import Job
from core.config import config

def count(
    input_fastqs,
    output,
    readset_name,
    ref_dir
    ):

    return Job(
        input_fastqs,
        [
            output 
        ],
        [
            ['cellranger_count', 'module_cellranger'] 
        ],
        command="""\
mkdir -p {folder} && \\
cd {folder} && \\
cellranger count \\
  --id={id} \\
  --transcriptome={ref} \\
  --fastqs={input} \\
  --sample={sample} \\
  {other_options}""".format(
            folder=os.path.dirname(os.path.dirname(os.path.dirname(output))),
            id=os.path.basename(os.path.dirname(os.path.dirname(output))),
            ref=ref_dir,
            input=os.path.dirname(input_fastqs[0]),
            sample=readset_name,
            other_options=config.param('cellranger_count', 'cellranger_other_options', required=False)
        )
    )

