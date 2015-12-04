#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def decompose_and_normalize_mnps(input, vt_output):

    return Job(
        [input],
        [vt_output],
        [
            ['DEFAULT', 'module_htslib'],
            ['decompose_and_normalize_mnps', 'module_vt']
        ],
        command="""\
zcat {input} | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r {reference_sequence} - | bgzip -cf > {vt_output} \\
        """.format(
        input=input,
        reference_sequence=config.param('DEFAULT', 'genome_fasta', type='filepath'),
        vt_output=vt_output
        )
    )
