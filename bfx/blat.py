#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def blat_dna_vs_dna(ref, query, output, other_options=""):
    return Job(
        [ref, query],
        [output],
        [['DEFAULT', 'module_ucsc']],
        command="""\
blat -t=dna -q=dna \\
  {ref} \\
  {query} \\
  {output} \\
  {other_options}""".format(
        ref=ref,
        query=query,
        output=output,
        other_options=other_options
        )
    )
