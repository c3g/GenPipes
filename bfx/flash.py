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
    output_dir,
    fastq_output,
    readset_name,
    log_output,
    hist_output,
    flash_stat_file
    ):

    # Paired end reads
    inputs = [input1, input2]
    outputs = [fastq_output, log_output, hist_output]

    flash_min, flash_max = None, None
    if flash_stat_file:                     # If flash_stat_file is present, then it means this is the Flash 2nd pass
        inputs.append(flash_stat_file)      # so flash_stat_file should be added to the job input list to well handle job dependencies
        # Then parse the stat file to build a stat dictionary
        flash_stats = dict()
        with open(flash_stat_file, 'r') as f:
            for line in f:
                splitLine = line.split()
                flash_stats[splitLine[1]] = ",".join(splitLine[2:3])
        # Finally retrieve from the dictionary the min & max overlap for the current readset
        [flash_min, flash_max] = flash_stats[readset_name].split(",")

    return Job(
        inputs,
        outputs,
        [
            ['flash', 'module_flash']
        ],
        command="""\
$FLASH_HOME/flash \\
  -t {threads} \\
  -m {min_overlap} \\
  -M {max_overlap} \\
  -o {name_out} \\
  {inputs} 2>&1 | tee {log_out}""".format(
        threads=config.param('flash', 'threads', type='posint'),
        min_overlap=flash_min if flash_min else config.param('flash', 'min_overlap', type='posint'),
        max_overlap=flash_max if flash_max else config.param('flash', 'max_overlap', type='posint'),
        name_out=os.path.join(output_dir, readset_name),
        inputs=" \\\n  ".join(inputs[0:2]),
        log_out=log_output
        ),
        removable_files=[fastq_output]
    )
