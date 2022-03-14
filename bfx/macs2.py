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

# MUGQIC Modules
from core.config import *
from core.job import *

def callpeak (options, genome_size, treatment_files, control_files, output_prefix_name, output, ini_section='macs2_callpeak'):

    return Job(
            treatment_files + control_files,
            output,
            [
                [ini_section, 'module_python'],
                [ini_section, 'module_macs2']
            ],
            command="""\
macs2 callpeak {options} {other_options} \\
  --tempdir {tmp_dir} \\
  --gsize {genome_size} \\
  --treatment \\
  {treatment_files}{control_files} \\
  --name {output_prefix_name} \\
  >& {output_prefix_name}.diag.macs.out""".format(
                options=options,
                other_options=global_config_parser.param(ini_section, 'other_options') if global_config_parser.param(ini_section, 'other_options') else "",
                tmp_dir=global_config_parser.param(ini_section, "tmp_dir"),
                genome_size=genome_size,
                treatment_files=" \\\n  ".join(treatment_files),
                control_files=" \\\n  --control \\\n  " + " \\\n  ".join(control_files) if control_files else " \\\n  --nolambda",
                output_prefix_name=output_prefix_name
            ),
            name="macs2_callpeak"
        )
