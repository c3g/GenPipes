#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def blat_dna_vs_dna(ref, query, outout, other_options=""):
    job = Job([ref, query], [outout], [['DEFAULT', 'module_ucsc']])

    job.command = \
"""blat -t=dna -q=dna \\
   {ref} \\
   {query} \\
   {output} \\
   {other_options} \\
   """.format(
        ref=ref,
        query=query,
        output=output,
        other_options=other_options
    )
    return job
