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
from core.config import *
from core.job import *


def differential_analysis(design_file, input_files, outputfiles, output_dir):

    suffix = re.sub(".*\.readset_", ".readset_", input_files[0])

    return Job(
        input_files,
        outputfiles,
        [
            ['methlykit_differential_analysis', 'module_R'],
            ['methlykit_differential_analysis', 'module_mugqic_tools']
        ],
        command="""\
Rscript $R_TOOLS/methylKit.R \\
  -design {design_file} \\
  -outdir {output_folder} \\
  -build {genome} \\
  -suff {input_suffix} \\
  {other_options}""".format(
            design_file=design_file,
            genome=config.param('methlykit_differential_analysis', 'assembly'),
            output_folder=output_dir,
            input_suffix=suffix,
            other_options=config.param('methlykit_differential_analysis', 'other_options')
        )
    )
    
