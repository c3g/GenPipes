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
    inputs,
    output,
    sample_id,
    sample_name,
    project,
    ref_dir,
    ini_section='cellranger_count'
    ):

    return Job(
        inputs,
        [output],
        [
            [ini_section, 'module_cellranger'] 
        ],
        command="""\
mkdir -p {folder} && \\
cd {folder} && \\
cellranger count \\
  --id={id} \\
  --transcriptome={ref} \\
  --fastqs={input} \\
  --project={project} \\
  --sample={sample} \\
  {other_options}""".format(
            folder=os.path.dirname(os.path.dirname(os.path.dirname(output))),
            id=sample_id,
            project=project,
            ref=ref_dir,
            input=os.path.dirname(inputs[0]),
            sample=sample_name,
            other_options=config.param(ini_section, 'other_options', required=False)
        )
    )

def atac(
    inputs,
    output,
    sample_id,
    sample_name,
    project,
    ref_dir,
    ini_section='cellranger_atac'
    ):

    return Job(
        inputs,
        [output],
        [
            [ini_section, 'module_cellranger_atac'] 
        ],
        command="""\
mkdir -p {folder} && \\
cd {folder} && \\
cellranger-atac count \\
  --id={id} \\
  --reference={ref} \\
  --fastqs={input} \\
  --project={project} \\
  --sample={sample} \\
  {other_options}""".format(
            folder=os.path.dirname(os.path.dirname(os.path.dirname(output))),
            id=sample_id,
            project=project,
            ref=ref_dir,
            input=os.path.dirname(inputs[0]),
            sample=sample_name,
            other_options=config.param(ini_section, 'other_options', required=False)
        )
    )

def vdj(
    inputs,
    output,
    sample_id,
    sample_name,
    project,
    ref_dir,
    ini_section='cellranger_vdj'
    ):

    return Job(
        inputs,
        [output],
        [
            [ini_section, 'module_cellranger'] 
        ],
        command="""\
mkdir -p {folder} && \\
cd {folder} && \\
cellranger vdj \\
  --id={id} \\
  --reference={ref} \\
  --fastqs={input} \\
  --project={project} \\
  --sample={sample} \\
  {other_options}""".format(
            folder=os.path.dirname(os.path.dirname(os.path.dirname(output))),
            id=sample_id,
            project=project,
            ref=ref_dir,
            input=os.path.dirname(inputs[0]),
            sample=sample_name,
            other_options=config.param(ini_section, 'other_options', required=False)
        )
    )