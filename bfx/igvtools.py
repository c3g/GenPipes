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

def compute_tdf(input, output):
    return Job(
        [input],
        [output],
        [
            ['compute_tdf', 'module_java'],
            ['compute_tdf', 'module_igvtools']
        ],
        command="""\
 java -Xmx{ram}  -Djava.awt.headless=true -jar $IGVTOOLS_JAR count {option} \\
  {input} \\
  {output} \\
  {genome}""".format(
        ram=config.param('igvtools_compute_tdf', 'ram'),
        option=config.param('igvtools_compute_tdf', 'option'),
        input=input,
        output=output,
        genome=config.param('compute_tdf', 'igv_genome', param_type='filepath')
        )
    )
