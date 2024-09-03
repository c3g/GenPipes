################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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
import logging
import os

from ..core.job import Job

log = logging.getLogger(__name__)

def run(input_file, pipeline, samples, readsets, job_name, metrics):
    """
    Calls job2json_project_tracking within jobs to update metrics
    """
    # The project tracking json file is provided as an environment variable.
    # Variable is exported prior to job submission to avoid having the json filename in the job script.
    # The json filename contains a timestamp that caused unwanted restarts.
    return Job(
        [input_file],
        [],
        [],
        command="""\
module purge && \\
{job2json_project_tracking_script} \\
  -s {samples} \\
  -r {readsets} \\
  -j {job_name} \\
  -o $PT_JSON_OUTFILE \\
  -m {metrics}""".format(
    job2json_project_tracking_script="job2json_project_tracking.py",
    samples=samples,
    readsets=readsets,
    job_name=job_name,
    metrics=metrics
    )
  )
