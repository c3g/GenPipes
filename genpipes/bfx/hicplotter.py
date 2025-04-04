################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
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
            hmc=global_conf.global_get('interaction_matrices_Chr', 'hmc'),
            fileNamePlot=fileNamePlot,
            newFileNamePlot=newFileNamePlot
        ),
        name=f"interaction_matrices_Chr.plotting.{sample_name}_{chrom}_res{res}"
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
            hmc=global_conf.global_get('interaction_matrices_Chr', 'hmc')
        ),
        name=f"interaction_matrices_genome.plotting.{sample_name}_res{res}"
    )
