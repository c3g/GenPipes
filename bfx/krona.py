#!/usr/bin/env python

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
