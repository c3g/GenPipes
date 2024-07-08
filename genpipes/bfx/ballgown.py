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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def ballgown( 
    input_files,
    design_file,
    output_dir,
    outputs
    ):

    return  Job(
        input_files + [design_file],
        outputs,
        [
            ['ballgown', 'module_mugqic_tools'],
            ['ballgown', 'module_R']
        ],
        command="""\
Rscript $R_TOOLS/ballgown.R \\
  -d {design_file} \\
  -o {output_dir} \\
  """.format(
            design_file=design_file,
            output_dir=output_dir,
        )
    )

