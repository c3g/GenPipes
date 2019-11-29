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

def calculate_reproducible_score( output_dir, sample1, sample2, file1_path, file2_path,  chromosome ):

    output_file= "".join((output_dir, "_".join(("/hicrep", sample1,  "vs", sample2 , chromosome)), ".tmp"))
    return Job(
        [file1_path,file2_path],
        [output_file],
        [
            ['reproducibility_scores', 'module_mugqic_tools'],
            ['reproducibility_scores', 'module_R']
        ],
        command="""\
        mkdir -p {output_dir} &&
Rscript /home/pubudu/projects/rrg-bourqueg-ad/pubudu/job_outputs/hicrep.R \\
  -s1 {sample1} \\
  -s2 {sample2} \\
  -f1 {file1_path} \\
  -f2 {file2_path} \\
  -c {chromosome} \\
  -o {output_dir}""".format(
        sample1=sample1,
        sample2=sample2,
        file1_path=file1_path,
        file2_path=file2_path,
        output_dir=output_dir,
        chromosome=chromosome
    ))


def merge_files( input_files, temp_out_dir):
    return Job(
        input_files,
        ["output_file"],
        [
            ['reproducibility_scores', 'module_mugqic_tools'],
            ['reproducibility_scores', 'module_R']
        ],
        command="""\
        cd {temp_out_dir} &&
        for d in */ ; do \\
        cat "$d"hicrep*.tmp | awk -v FS='\\t' -v OFS='\\t' -v dir="${{d%/*}}" '{{sum+=$1}} END {{print dir,sum/NR}}' > "${{d%/*}}".tmp \\
        done &&
        cat *.tmp > hicrep_reproducibilityscore_matrix.tsv &&
        rm *.tmp""".format(
        temp_out_dir=temp_out_dir
        ))