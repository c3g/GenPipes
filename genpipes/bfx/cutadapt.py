################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

log = logging.getLogger(__name__)

def trim(input1, input2, prefix, adapter_3p_fwd, adapter_3p_rev):

    if input2:  # Paired end reads
        output_pair1 = prefix + ".trim.pair1.fastq.gz"
        output_pair2 = prefix + ".trim.pair2.fastq.gz"
        inputs = [input1, input2]
        output = [output_pair1, output_pair2]
    else:   # Single end reads
        output_pair1 = prefix + ".trim.single.fastq.gz"
        inputs = [input1]
        output = [output_pair1]

    return Job(
        inputs,
        output,
        [
            ['cutadapt', 'module_cutadapt']
        ],

        command="""\
cutadapt {adapter_5p_fwd} \\
  {adapter_5p_rev} \\
  {adapter_3p_fwd} \\
  {adapter_3p_rev} \\
  {nthread} \\
  {options} \\
  {output_fwd} \\
  {output_rev} \\
  {inputs}""".format(
      adapter_5p_fwd="-g " + global_conf.global_get('cutadapt', 'adapter_5p_fwd') if global_conf.global_get('cutadapt', 'adapter_5p_fwd') else "",
      adapter_5p_rev="-G " + global_conf.global_get('cutadapt', 'adapter_5p_rev') if input2 and global_conf.global_get('cutadapt', 'adapter_5p_fwd') else "",
      adapter_3p_fwd="-a " + adapter_3p_fwd,
      adapter_3p_rev="-A " + adapter_3p_rev if input2 else "",
      options=global_conf.global_get('cutadapt', 'options'),
      nthread="-j " + global_conf.global_get('cutadapt', 'threads'),
      inputs=" \\\n  ".join(inputs),
      output_fwd="-o " + output_pair1,
      output_rev="-p " + output_pair2 if input2 else ""
      ),
    )
