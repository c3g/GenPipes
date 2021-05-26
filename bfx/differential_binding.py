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

def diffbind2(input_files, comparison, design, output_file):

    #merge all the tsvs and create final .csv file

    if not isinstance(input_files, list):
       input_files = [input_files]

    return Job(
        input_files,
        [output_file],
        [
            ['macs2_callpeak', 'module_macs2']
        ],
        command="""\
                    mkdir -p {output_dir}/""".format(
            output_dir=output_file
        ))
def diffbind( input_files, comparison, design, readset, output_dir):

    output_file =  "".join((output_dir, "_".join(("/diffbind",comparison,"dba.txt"))))
    return Job(
        input_files,
        [output_file],
        [
            ['differential_binding', 'module_mugqic_tools'],
            ['differential_binding', 'module_R']
        ],
        command="""\
        mkdir -p {output_dir} &&
#Rscript $R_TOOLS/diffbind.R \\
Rscript /home/pubudu/projects/rrg-bourqueg-ad/pubudu/chipseq_diff/diffbind.R \\
  -d {design} \\
  -r {readset} \\  
  -c {comparison} \\
  -o {output_file}""".format(
        design=design,
        comparison=comparison,
        output_file=output_file,
        output_dir=output_dir,
        readset=readset
    ))
