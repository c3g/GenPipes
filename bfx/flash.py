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
        threads=global_config_parser.param('flash', 'threads', param_type='posint'),
        min_overlap="${minFlashOverlap}" if pre_command else global_config_parser.param('flash', 'min_overlap', param_type='posint'),
        max_overlap="${maxFlashOverlap}" if pre_command else global_config_parser.param('flash', 'max_overlap', param_type='posint'),
        name_out=re.sub(".extendedFrags.fastq", "", fastq_output),
        inputs=" \\\n  ".join(inputs[0:2]),
        log_out=log_output
        ),
        removable_files=[fastq_output]
    )
