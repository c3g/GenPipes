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

from core.job import *

log = logging.getLogger(__name__)

def run(input_file, file_regex=None):
    """
    Calls job2json_project_tracking_rm to remove deleted files when cleaning up GenPipes output files.
    """

    if not file_regex:
        file_regex = input_file
    return Job(
        [input_file],
        [],
        [],
        command="""\
module load {module_python} && \\
{job2json_project_tracking_rm_script} \\
  -f {file_regex} \\
  -o $PT_JSON_OUTFILE \\
&& \\
module unload {module_python}""".format(
    module_python=config.param('DEFAULT', 'module_python'),
    job2json_project_tracking_rm_script=os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "utils", "job2json_project_tracking_rm.py"),
    file_regex=file_regex
    )
  )
