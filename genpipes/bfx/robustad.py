################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
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

def call_TADs(matrix, output_dir, res, minwin, maxwin):

    prefix = os.path.splitext(os.path.basename(matrix))[0]

    output_Scores = os.path.join(output_dir, "".join(("BoundaryScores_", prefix, "_binSize" , f'{float(int(res)/1000):g}' ,"_minW" + minwin + "_maxW" + maxwin + "_minRatio1.5.txt")))
    output_calls = os.path.join(output_dir, "".join(("TADBoundaryCalls_", prefix, "_binSize" , f'{float(int(res)/1000):g}' ,"_minW" + minwin + "_maxW" + maxwin + "_minRatio1.5_threshold0.2.txt")))

    return Job(
        [matrix],
        [output_Scores, output_calls],
        [
            ["identify_TADs", "module_R"],
            ["identify_TADs", "module_mugqic_tools"]
        ],
        command="""
Rscript {RobusTAD} \\
  -i {input_matrix} \\
  -H \\
  -b {res} \\
  --minWin {minwin} \\
  --maxWin {maxwin} \\
  -o {ouput_dir}""".format(
            ouput_dir=output_dir,
            RobusTAD=os.path.expandvars("${R_TOOLS}/RobusTAD.R"),
            input_matrix=matrix,
            res=int(res)/1000,
            minwin=minwin,
            maxwin=maxwin
        ),
        removable_files=[matrix]
    )
