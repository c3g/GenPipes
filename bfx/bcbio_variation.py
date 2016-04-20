#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
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
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def tumor_pair_ensemble(input_mutect2, input_vardict, input_samtools, output):
    
    return Job(
        [input_mutect2,input_vardict,input_samtools],
        [output],
        [
            ['bcbio_ensemble', 'module_java'],
            ['bcbio_ensemble', 'module_bcbio_variation']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $BCBIO_VARIATION_JAR \\
  $BCBIO_VARIATION_HOME/tumor_pair_ensemble.yaml \\
  {reference_sequence} \\
  {output} \\
  {input_mutect2} {input_vardict} {input_samtools}""".format(
        tmp_dir=config.param('bcbio_ensemble', 'tmp_dir'),
        java_other_options=config.param('bcbio_ensemble', 'java_other_options'),
        ram=config.param('bcbio_ensemble', 'ram'),
        reference_sequence=config.param('bcbio_ensemble', 'genome_fasta', type='filepath'),
        output=output,
        input_mutect2=input_mutect2,
        input_vardict=input_vardict,
        input_samtools=input_samtools
        )
    )

