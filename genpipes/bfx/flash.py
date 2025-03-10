################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import logging
import re

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

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
    min_overlap = global_conf.global_get('flash', 'min_overlap', param_type='posint')
    max_overlap = global_conf.global_get('flash', 'max_overlap', param_type='posint')
    if flash_stats_file:                        # If flash_stat_file is present, then it means this is the Flash 2nd pass
        inputs.append(flash_stats_file)         # Append it to inputs for dependencies matter
        pre_command = f"""\
minFlashOverlap=$(grep -w {readset_name} {flash_stats_file} | cut -f 5)
maxFlashOverlap=$(grep -w {readset_name} {flash_stats_file} | cut -f 6)"""
        min_overlap = "${minFlashOverlap}"
        max_overlap = "${maxFlashOverlap}"

    return Job(
        inputs,
        outputs,
        [
            ['flash', 'module_flash']
        ],
        command=f"""\
{pre_command}
$FLASH_HOME/flash \\
  -t {global_conf.global_get('flash', 'threads', param_type='posint')} \\
  -m {min_overlap} \\
  -M {max_overlap} \\
  -o {re.sub(".extendedFrags.fastq", "", fastq_output)} \\
  {" \\\n ".join(inputs[0:2])} 2>&1 | tee {log_output}""",
        removable_files=[fastq_output]
    )
