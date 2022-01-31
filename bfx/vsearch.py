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
import logging

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

def uchime(
    cat_sequence_fasta,
    filter_fasta
    ):

    inputs = [cat_sequence_fasta]
    outputs = [filter_fasta]

    return Job(
        inputs,
        outputs,
        [
            ['uchime', 'module_vsearch'],
        ],
        command="""\
  $VSEARCH_HOME/usearch61 \\
  --uchime_ref {cat_sequence_fasta} \\
  --db {database} \\
  --nonchimeras {filter_fasta} \\
  --threads {threads_number}""".format(
        cat_sequence_fasta=cat_sequence_fasta,
        database=config.param('uchime', 'chimera_database'),
        threads_number=config.param('uchime', 'threads'),
        filter_fasta=filter_fasta
        ),
        removable_files=[filter_fasta]
    )
