################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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
import os

log = logging.getLogger(__name__)

def run(pipeline, samples, readsets, job_name, metrics):
    """
    Calls job2json_project_tracking within jobs to update metrics
    """
    pipeline_output_dir = pipeline.output_dir
    json_folder = os.path.join(pipeline_output_dir, "json")
    operation_name = pipeline.__class__.__name__
    timestamp = pipeline.timestamp
    json_outfile = os.path.join(json_folder, f"{operation_name}_{timestamp}.json")

    return """\
{job2json_project_tracking_script} \\
  -s \\"{samples}\\" \\
  -r \\"{readsets}\\" \\
  -j \\"{job_name}\\" \\
  -o \\"{json_outfile}\\" \\
  -m {metrics}""".format(
    job2json_project_tracking_script=os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "utils", "job2json_project_tracking.py"),
    samples=",".join([sample.name for sample in samples]),
    readsets=",".join([readset.name for readset in readsets]),
    job_name=job_name,
    metrics=metrics,
    json_outfile=json_outfile
  )
