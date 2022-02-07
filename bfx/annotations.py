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

def structural_variants(input, output):
    return Job(
        [input],
        [os.path.abspath(output)],
        [
            ['sv_annotation', 'module_python'],
            ['sv_annotation', 'module_sv_annotations']
        ],
        command="""\
rm -f {output} && \\
python $SVANNOT_PATH/simple_sv_annotation.py \\
        -g $SVANNOT_PATH/az-cancer-panel.txt \\
        -k $SVANNOT_PATH/fusion_pairs.txt  \\
        -o {output} \\
        {input}""".format(
            input=input,
            output=output,
        )
 )
