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
from core.job import *

log = logging.getLogger(__name__)

def krona(
    otu_normalized_table,
    sample_name,
    alpha_diversity_krona_file
    ):

    inputs = [otu_normalized_table]
    outputs = [alpha_diversity_krona_file]

    return Job(
        inputs,
        outputs,
        [
            ['krona', 'module_perl'],
            ['krona', 'module_python'],
            ['krona', 'module_krona']
        ],

        command="""\
cat {sample_name} > _temp && \\
$PERL_HOME/bin/perl $KRONA_HOME/ImportText.pl \\
  _temp \\
  -o {alpha_diversity_krona_file}
rm -f _temp""".format(
        sample_name=' '.join(sample_name),
        alpha_diversity_krona_file=alpha_diversity_krona_file
        ),
        removable_files=[alpha_diversity_krona_file]
    )
