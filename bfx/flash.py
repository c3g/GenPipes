#!/usr/bin/env python

# Python Standard Modules
import logging
import os
import gzip

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

#the merge file is name.extendedFrags.fastq 
def flash(
    input1,
    input2,
    fastq_output,
    readset_name,
    log_output,
    hist_output,
    flash_stats_file
    ):

    # Paired end reads
    inputs = [input1, input2]
    outputs = [fastq_output, log_output, hist_output]

    pre_command = None
    if flash_stats_file:                        # If flash_stat_file is present, then it means this is the Flash 2nd pass
        inputs.append(flash_stats_file)         # Append it to inputs for dependencies matter
        pre_command="""\
minFlashOverlap=$(grep {readset} {file} | cut -f 5)
maxFlashOverlap=$(grep {readset} {file} | cut -f 6)""".format(
            readset=readset_name,
            file=flash_stats_file
        )

    return Job(
        inputs,
        outputs,
        [
            ['flash', 'module_flash']
        ],
        command="""\
{pre_command}
$FLASH_HOME/flash \\
  -t {threads} \\
  -m {min_overlap} \\
  -M {max_overlap} \\
  -o {name_out} \\
  {inputs} 2>&1 | tee {log_out}""".format(
        pre_command=pre_command,
        threads=config.param('flash', 'threads', type='posint'),
        min_overlap="${minFlashOverlap}" if pre_command else config.param('flash', 'min_overlap', type='posint'),
        max_overlap="${maxFlashOverlap}" if pre_command else config.param('flash', 'max_overlap', type='posint'),
        name_out=re.sub(".extendedFrags.fastq", "", fastq_output),
        inputs=" \\\n  ".join(inputs[0:2]),
        log_out=log_output
        ),
        removable_files=[fastq_output]
    )
