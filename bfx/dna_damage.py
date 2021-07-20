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

def estimate_damage(mpileup_r1, mpileup_r2, output, sample_id, options):

    return Job(
        [mpileup_r1,mpileup_r2],
        [output],
        [
            ['DEFAULT', 'module_dna_damage'],
        ],
        command="""\
perl $DNA_DAMAGE_PATH/estimate_damage.pl {options} \\
    --mpileup1 {mpileup_r1} \\
    --mpileup2 {mpileup_r2} \\
    --id {sample_id} \\
    {output}""".format(
        options=options,
        mpileup_r1=mpileup_r1,
        mpileup_r2=mpileup_r2,
        sample_id=sample_id,
        output="> " + output if output else "",
        )
    )

def estimate_damage_r(input, output):

    return Job(
        [input],
        [output],
        [
            ['DEFAULT', 'module_dna_damage'],
            ['DEFAULT', 'module_R'],
        ],
        command="""\
Rscript --vanilla $DNA_DAMAGE_PATH/plot_damage.R \\
    {input} \\
    {output}""".format(
        input=input,
        output=output,
        )
    )
