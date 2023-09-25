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
import logging

# MUGQIC Modules
from core.config import *
from core.job import *

### Start from here ###

def paired(input1, input2, output_directory, sample):
	return Job(
		input_files=[input1, input2],
		output_files=[output_directory, sample],
		module_entries=[
            ['sortmerna', 'module_sortmerna']
        ],
        command="""\
sortmerna --threads {cpu} \\
	--ref $SORTMERNA_DATA/rRNA_databases/rfam-5.8s-database-id98.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/rfam-5s-database-id98.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/silva-arc-16s-id95.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/silva-arc-23s-id98.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/silva-bac-16s-id90.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/silva-bac-23s-id98.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/silva-euk-18s-id95.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/silva-euk-28s-id98.fasta \\
	--reads {input1} \\
	--reads {input2} \\
    --aligned {output_directory}/${sample}/${sample}.aligned \\
    --kvdb {output_directory}/${sample}/kvdb \\
    --readb {output_directory}/${sample}/readb \\
    --idx-dir {output_directory}/idx-dir""".format(
			input1=input1,
			input2=input2,
        	cpu=config.param('sortmerna', 'cluster_cpu', required=True), 
        	output_directory=output_directory,
        	sample=sample,
        	other_options=(" \\\n  " + config.param('sortmerna', 'other_options', required=False)) if config.param('sortmerna', 'other_options', required=False) else ""
        )
	)

def single(input1, output_directory, sample):
	return Job(
		input_files=[input1],
		output_files=[output_directory, sample],
		module_entries=[
            ['sortmerna', 'module_sortmerna'],
        ],
        command="""\
sortmerna --threads {cpu} \\
	--ref $SORTMERNA_DATA/rRNA_databases/rfam-5.8s-database-id98.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/rfam-5s-database-id98.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/silva-arc-16s-id95.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/silva-arc-23s-id98.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/silva-bac-16s-id90.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/silva-bac-23s-id98.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/silva-euk-18s-id95.fasta \\
	--ref $SORTMERNA_DATA/rRNA_databases/silva-euk-28s-id98.fasta \\
	--reads {input1} \\ 
    --aligned {output_directory}/${sample}/${sample}.aligned \\
    --kvdb {output_directory}/${sample}/kvdb \\
    --readb {output_directory}/${sample}/readb \\
    --idx-dir {output_directory}/idx-dir""".format(
			input1=input1,
        	cpu=config.param('sortmerna', 'cluster_cpu', required=True), 
        	output_directory=output_directory,
        	sample=sample,
        	other_options=(" \\\n  " + config.param('sortmerna', 'other_options', required=False)) if config.param('sortmerna', 'other_options', required=False) else ""
        )
	)