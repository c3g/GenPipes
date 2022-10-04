################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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

def calculate_reproducible_score(output_file, sample1, sample2, file1_path, file2_path, chromosome , resolution, bound_width, weights, corr, down_sampling, smooth):

    # output_file= "".join((output_dir, "_".join(("/hicrep", sample1,  "vs", sample2 , chromosome, resolution, "res",smooth,
                                                # bound_width, down_sampling)), ".tmp"))

    return Job(
        [file1_path, file2_path],
        [output_file],
        [
            ['reproducibility_scores', 'module_mugqic_tools'],
            ['reproducibility_scores', 'module_R']
        ],
        command="""\
Rscript $R_TOOLS/hicrep.R \\
  -s1 {sample1} \\
  -s2 {sample2} \\
  -f1 {file1_path} \\
  -f2 {file2_path} \\
  -c {chromosome} \\
  -o {output_file} \\
  -r {resolution} \\
  -sm {smooth} \\
  -b {bound_width} \\
  -w {weights} \\
  -cor {corr} \\
  -d {down_sampling}""".format(
        sample1=sample1,
        sample2=sample2,
        file1_path=file1_path,
        file2_path=file2_path,
        output_file=output_file,
        chromosome=chromosome,
        resolution=resolution,
        smooth=smooth,
        bound_width=bound_width,
        weights=weights,
        corr=corr,
        down_sampling=down_sampling
        )
  )


def merge_tmp_files(input_files, output_files, temp_out_dir, resolution, smooth, bound_width, down_sampling, temp_dir):
    if not isinstance(input_files, list):
        input_files = [input_files]
    if not isinstance(output_files, list):
        output_files = [output_files]

    return Job(
        input_files,
        output_files,
        [
            ['reproducibility_scores', 'module_mugqic_tools']
        ],
        command="""\\
cd {temp_out_dir}/{temp_dir} &&
for d in */ ; do 
    awk -v OFS="_" 'NR==0 {{print; next}} FNR==0 {{next}} {{print FILENAME, $0}}' "$d"hicrep*_{resolution}_res_{smooth}_{bound_width}_{down_sampling}*.tmp |  \\
    awk -F "_" 'BEGIN{{printf "chr\\tSCC\\tSD\\tSmoothing\\n"}}{{printf $(NF-6)"\\t"$NF"\\n"}}' > "${{d%/*}}"_res_{resolution}_{smooth}_{bound_width}_{down_sampling}.tsv 
done""".format(
                temp_out_dir=temp_out_dir,
                resolution=resolution,
                smooth=smooth,
                bound_width=bound_width,
                down_sampling=down_sampling,
                temp_dir=temp_dir
                )
            )

def merge_tsv(input_files, out_dir, output_file, temp_dir):

    #merge all the tsvs and create final .csv file

    if not isinstance(input_files, list):
        input_files = [input_files]

    return Job(
        input_files,
        [os.path.join(out_dir, output_file)],
        [
            ['reproducibility_scores', 'module_mugqic_tools']
        ],
        command="""\\
cd {out_dir}/{temp_dir} && \\
awk -v OFS="_" 'NR==1 {{print; next}} FNR==1 {{next}} {{print substr( FILENAME, 1, length(FILENAME)-4 ),$0}}' *.tsv | tr "\\t" "," > ../{output_file}""".format(
                        out_dir=out_dir,
                        output_file=output_file,
                        temp_dir=temp_dir
                        )
            )
