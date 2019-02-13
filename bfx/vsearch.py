#!/usr/bin/env python

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
