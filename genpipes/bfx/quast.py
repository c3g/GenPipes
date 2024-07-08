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
import os
import logging

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

log = logging.getLogger(__name__)

def quast(input, output_dir, ini_section='quast_consensus_metrics'):
    inputs = [input]
    outputs = [
        os.path.join(output_dir, "report.html"),
        os.path.join(output_dir, "report.pdf"),
        os.path.join(output_dir, "report.tex"),
        os.path.join(output_dir, "report.tsv"),
        os.path.join(output_dir, "report.txt")
        ]

    return Job(
        inputs,
        outputs,
        [
            ['quast', 'module_quast']
        ],

        command="""\
export LC_ALL=en_CA.UTF-8 && \\
quast.py {reference} \\
  {output_dir} \\
  {features} \\
  {nthread} \\
  {input}""".format(
      reference="-r " + global_conf.global_get(ini_section, 'reference_genome') if global_conf.global_get(ini_section, 'reference_genome', required=False) else "",
      features="--features " + global_conf.global_get(ini_section, 'genomic_feature', required=False) if global_conf.global_get(ini_section, 'genomic_feature', required=False) else "",
      nthread="--threads " + global_conf.global_get(ini_section, 'threads') if global_conf.global_get(ini_section, 'threads', required=False) else "",
      output_dir="--output-dir " + output_dir,
      input=input
      ),
    )
