################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def run(input1, input2, output_directory, output_directory_sample, sample, ini_section="sortmerna"):
	return Job(
		input_files=[input1, input2],
		output_files=[output_directory_sample, os.path.join(output_directory_sample, sample + ".aligned.log")],
		module_entries=[
            [ini_section, 'module_sortmerna']
        ],
        command="""\
sortmerna --threads {threads} {other_options} \\
	{ref_rfam_58s} \\
	{ref_rfam_5s} \\
	{ref_silva_arc_16s} \\
	{ref_silva_arc_23s} \\
	{ref_silva_bac_16s} \\
	{ref_silva_bac_23s} \\
	{ref_silva_euk_18s} \\
	{ref_silva_euk_28s} \\
	--reads {input1} \\
	{input2} \\
    --aligned {output_directory_sample}/{sample}.aligned \\
    --kvdb {output_directory_sample}/kvdb \\
    --readb {output_directory_sample}/readb""".format(
			input1=input1,
			input2="--reads " + input2 if input2 else "",
        	ref_rfam_58s="--ref " + config.param(ini_section, 'rfam_5.8s', required=False) if config.param(ini_section, 'rfam_5.8s', required=False) else "",
			ref_rfam_5s="--ref " + config.param(ini_section, 'rfam_5s', required=False) if config.param(ini_section, 'rfam_5s', required=False) else "",
			ref_silva_arc_16s="--ref " + config.param(ini_section, 'silva_arc_16s', required=False) if config.param(ini_section, 'silva_arc_16s', required=False) else "",
			ref_silva_arc_23s="--ref " + config.param(ini_section, 'silva_arc_23s', required=False) if config.param(ini_section, 'silva_arc_23s', required=False) else "",
			ref_silva_bac_16s="--ref " + config.param(ini_section, 'silva_bac_16s', required=False) if config.param(ini_section, 'silva_bac_16s', required=False) else "",
			ref_silva_bac_23s="--ref " + config.param(ini_section, 'silva_bac_23s', required=False) if config.param(ini_section, 'silva_bac_23s', required=False) else "",
			ref_silva_euk_18s="--ref " + config.param(ini_section, 'silva_euk_18s', required=False) if config.param(ini_section, 'silva_euk_18s', required=False) else "",
			ref_silva_euk_28s="--ref " + config.param(ini_section, 'silva_euk_28s', required=False) if config.param(ini_section, 'silva_euk_28s', required=False) else "",
			threads=config.param(ini_section, 'threads', required=True), 
        	output_directory=output_directory,
        	output_directory_sample=output_directory_sample,
        	sample=sample,
        	other_options = config.param(ini_section, 'other_options')
        )
	)