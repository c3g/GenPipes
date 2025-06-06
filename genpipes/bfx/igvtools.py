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
        ram=global_conf.global_get('igvtools_compute_tdf', 'ram'),
        option=global_conf.global_get('igvtools_compute_tdf', 'option'),
        input=input,
        output=output,
        genome=global_conf.global_get('compute_tdf', 'igv_genome', param_type='filepath')
        )
    )
