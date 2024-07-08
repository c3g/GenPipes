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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

log = logging.getLogger(__name__)

def run(input_folder, output_folder, sample_name):
    input = os.path.join(input_folder, sample_name + ".arriba.tsv")
    #output = os.path.join(output_folder, sample_name + ".metacaller.1.tsv")
    output = os.path.join(output_folder, sample_name + ".metacaller.2.tsv")
    return Job(
        [input],
        [output],
        [
            ['fusionmetacaller', 'module_mugqic_tools'],
            ['fusionmetacaller', 'module_R'],
        ],
        command="""\\
Rscript $R_TOOLS/RunFusionMetaCaller.R \\
    {INPUT_FOLDER}   \\
    {OUTPUT_BASE_NAME}""".format(
            INPUT_FOLDER=input_folder,
            OUTPUT_BASE_NAME=os.path.join(output_folder, sample_name)
        )
    )
