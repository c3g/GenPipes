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

# MUGQIC Modules
from core.config import *
from core.job import *

def run(normal, tumor, normal_name, tumor_name, output_dir, other_options=None):
    normal_output = os.path.join(output_dir, normal_name + ".cobalt.ratio.pcf")
    tumor_output = os.path.join(output_dir, tumor_name + ".cobalt.ratio.pcf")

    return Job(
        [normal, tumor],
        [normal_output, tumor_output],
        [
            ['cobalt', 'module_java'],
            ['cobalt', 'module_R'],
            ['cobalt', 'module_cobalt'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $COBALT_JAR \\
  -threads {threads} \\
  -reference {reference} \\
  -reference_bam {reference_bam} \\
  -tumor {tumor} \\
  -tumor_bam {tumor_bam} \\
  -gc_profile {gc_profile} \\
  -output_dir {output_dir}""".format(
        tmp_dir=global_config_parser.param('cobalt', 'tmp_dir'),
        java_other_options=global_config_parser.param('cobalt', 'java_other_options'),
        ram=global_config_parser.param('cobalt', 'ram'),
        threads=global_config_parser.param('cobalt', 'threads'),
        gc_profile=global_config_parser.param('cobalt', 'gc_profile'),
        reference=normal_name,
        reference_bam=normal,
        tumor=tumor_name,
        tumor_bam=tumor,
        other_options=other_options,
        output_dir=output_dir,
        )
    )
