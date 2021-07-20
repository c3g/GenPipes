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

from core.config import *
from core.job import *

#Assuming that the module file has already been written and the scripts have been added to path

def run(sample_name, normal_name, tumor_name):
    normal_inputs = [os.path.join("pairedVariants", sample_name, "sCNAphase", normal_name + ".chr" + str(chr) + ".haps") for chr in range(1, 22)]
    #tumor_inputs = [os.path.join(tumor_name + ".chr" + str(chr) + ".vcf") for chr in range(1, 22)]

    return Job(
        normal_inputs,
        [None],
        [
            ['scnaphase', 'module_mugqic_tools'],
            ['scnaphase', 'module_R'],
        ],
        command="""\\
Rscript $R_TOOLS/RunsCNAphaseAnalysis.R \\
    {output_base_name} \\
    {normal_name}   \\
    {tumor_name}""".format(
            normal_name=normal_name,
            tumor_name=tumor_name,
            output_base_name=sample_name + ".inferCN",
        )
    )
