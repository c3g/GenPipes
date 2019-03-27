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

#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def decompose_and_normalize_mnps(inputs, vt_output=None):
    if not isinstance(inputs, list):
        inputs = [inputs]

    return Job(
        inputs,
        [vt_output],
        [
            ['decompose_and_normalize_mnps', 'module_htslib'],
            ['decompose_and_normalize_mnps', 'module_vt']
        ],
        command="""\
zless {input} | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r {reference_sequence} -  \\
        {vt_output}""".format(
        input=" \\\n  ".join(input for input in inputs),
        reference_sequence=config.param('decompose_and_normalize_mnps', 'genome_fasta', type='filepath'),
        vt_output="> " + vt_output if vt_output else " ",
        )
    )


def decompose_uniq_and_normalize_mnps(input_vcf, vt_output):
        ## iupac code is not supported by bcbio_variation_recall, remove variants with these code
        iupac_ambigous_code = ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"]
        return Job(
                [input_vcf],
                [vt_output],
                [
                        ['decompose_and_normalize_mnps', 'module_vt']
                ],
                command="""\
cat {input_vcf} | sed 's/ID=AD,Number=./ID=AD,Number=R/' | awk '$1~/#/ || ($4$5)!~/{iupac_str}/' | vt decompose -s - | vt uniq - | vt normalize - -r {reference_sequence} -o {vt_output} \\
                """.format(
                input_vcf=input_vcf,
                reference_sequence=config.param('vt', 'genome_fasta', type='filepath'),
                vt_output=vt_output,
                iupac_str="|".join(iupac_ambigous_code)
                )
        )

def sort(input, output, options):

    return Job(
        [input],
        [output],
        [
            ['vt_sort', 'module_htslib'],
            ['vt_sort', 'module_vt']
        ],
        command="""\
vt sort {options} -o {output} {input}""".format(
        options=options,
        input=" \\\n " + input if input else "-",
        output=output
        )
    )


