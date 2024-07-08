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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def run(input_bam, input_gtf, output_directory):
    input_base = os.path.basename(input_bam)
    return Job(
        [input_bam, input_gtf],
        [os.path.join(output_directory, input_base + ".metrics.tsv")],
        [
            ['rnaseqc2', 'module_rnaseqc2']
        ],
        command="""\
rnaseqc \\
  {input_gtf} \\
  {input_bam} \\
  {other_options} \\
  {output_directory}""".format(
            input_gtf=input_gtf,
            input_bam=input_bam,
            output_directory=output_directory,
            other_options=global_conf.global_get('rnaseqc2', 'other_options', required=False)
        )
    )

def parse_expression_profiling_efficiency_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export expression_profiling_efficiency=`awk '{{if ($0 ~ /Expression Profiling Efficiency/) {{print $4}}}}' {input_file}`"""
        )

def parse_rrna_rate_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export rrna_rate=`awk '{{if ($0 ~ /rRNA Rate/) {{print $3}}}}' {input_file}`"""
        )
