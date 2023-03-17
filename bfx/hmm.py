################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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

#MUGQIC Modules
from core.config import *
from core.job import *

def readCounter(
        input,
        output
        ):
    return Job(
            [input],
            [output],
            [
                ['hmm_readCounter', 'module_hmm']
            ],

            command="""\
readCounter \\
    --window {window_size} \\
    --quality {threshold} \\
    --chromosome {chr_list} \\
    {input} \\
    > {output}""".format(
        window_size=config.param('hmm_readCounter', 'window_size'),
        threshold=config.param('hmm_readCounter', 'threshold'),
        chr_list=config.param('hmm_readCounter', 'chr_list'),
        input=input,
        output=output
        )
    )
