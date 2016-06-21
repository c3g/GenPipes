#!/usr/bin/env python

# Python Standard Modules
import logging
import os

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

#the merge file is name.extendedFrags.fastq 
def flash(
    input1,
    input2,
    dir_output,
    merge_output,
    name_output,
    log_output
    ):

    # Paired end reads
    inputs = [input1, input2]
    outputs = [merge_output, log_output]

    return Job(
        inputs,
        outputs,
        [
            ['flash', 'module_flash']
        ],

        # CAUTION: Trimmomatic settings order is IMPORTANT!
        # FIRST Illuminaclip settings, THEN headcrop length, THEN trailing min quality, THEN minimum length
        command="""\
  $FLASH_HOME/flash \\
  -t {threads} \\
  -m {min_overlap} \\
  -M {max_overlap} \\
  -o {name_out} \\
  {inputs} 2>&1 | tee {log_out}""".format(
        threads=config.param('flash', 'threads', type='posint'),
        min_overlap=config.param('flash', 'min_overlap', type='posint'),
        max_overlap=config.param('flash', 'max_overlap', type='posint'),
        name_out=os.path.join(dir_output,name_output),
        inputs=" \\\n  ".join(inputs),
        log_out=log_output
        ),
        removable_files=[merge_output]
    )
