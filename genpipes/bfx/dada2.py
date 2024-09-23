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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def dada2(
    inputs,
    ampliconLengthFile,
    rawReadsFolder,
    designFile,
    output_directory
    ):
    inputs.append(ampliconLengthFile)
    return  Job(
        inputs,
        [output_directory],
        [
            ['dada2','module_R'],
            ['dada2','module_mugqic_tools'],
        ],
        command="""\
Rscript $R_TOOLS/asva.R \\
  -r {rawReadsFolder} \\
  -d {designFile} \\
  -o {output_directory} \\
  -tr {trainset} \\
  -tax {taxonomy} \\
  -p {pool_parameter} \\
  -amp {amplicon_length_file}""".format(
            rawReadsFolder=rawReadsFolder,
            designFile=designFile,
            output_directory=output_directory,
            trainset=global_conf.global_get('database', 'dada2_trainset'),
            taxonomy=global_conf.global_get('database', 'dada2_taxonomy'),
            pool_parameter=global_conf.global_get('dada2', 'pool_parameter'),
            amplicon_length_file=ampliconLengthFile
        )
    )

