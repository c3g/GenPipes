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

# MUGQIC Modules
from core.config import *
from core.job import *

def intra_chrom_matrix_plot(fileNameRN, sample_name, chrom, res, output_dir, fileNamePlot, newFileNamePlot):
    return Job(
        [fileNameRN],
        [newFileNamePlot],
        [
            ["interaction_matrices_Chr", "module_python"],
            ["interaction_matrices_Chr", "module_HiCPlotter"]
        ],
        command="""\
HiCPlotter.py \\
  -f {fileNameRN} \\
  -n {name} \\
  -chr {chr} \\
  -r {res} \\
  -fh 0 \\
  -o {output_dir} \\
  -ptr 0 \\
  -hmc {hmc} && \\
mv {fileNamePlot} {newFileNamePlot}""".format(
            res=res,
            chr=chrom,
            fileNameRN=fileNameRN,
            name=sample_name,
            output_dir=output_dir,
            hmc=global_config_parser.param('interaction_matrices_Chr', 'hmc'),
            fileNamePlot=fileNamePlot,
            newFileNamePlot=newFileNamePlot
        ),
        name="interaction_matrices_Chr.plotting."+sample_name+"_"+chrom+"_res"+res
    )

def genome_wide_matrix_plot(fileNameRN, sample_name, res, output_dir, output_file):
    return Job(
        [fileNameRN],
        [output_file],
        [
            ["interaction_matrices_Chr", "module_python"],
            ["interaction_matrices_Chr", "module_HiCPlotter"]
        ],
        command="""\
HiCPlotter.py \\
  -f {fileNameRN} \\
  -n {name} \\
  -chr Genome \\
  -r {res} \\
  -fh 0 \\
  -o {output_dir} \\
  -ptr 0 \\
  -hmc {hmc} \\
  -wg 1""".format(
            res=res,
            fileNameRN=fileNameRN,
            name=sample_name, 
            output_dir=output_dir,
            hmc=global_config_parser.param('interaction_matrices_Chr', 'hmc')
        ),
        name="interaction_matrices_genome.plotting."+sample_name+"_res"+res
    )
