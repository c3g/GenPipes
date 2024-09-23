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
import os

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def run(input, ped, output_dir):
	output = os.path.join(output_dir, ".ped_check.csv")
	return Job(
        [input, ped],
        [output],
        [
            ['run_peddy', 'module_python'],
        ],
        command="""\
python -m peddy --plot {options} \\
    --procs {threads} \\
    --prefix {output} {input} {ped}""".format(
            options=global_conf.global_get('run_peddy', 'options'),
	        threads=global_conf.global_get('run_peddy', 'threads'),
	        input=input,
	        ped=ped,
	        output=output_dir,
        )
    )
