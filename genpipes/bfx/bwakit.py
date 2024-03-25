################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
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
# along with GenPipes. If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import os

# GenPipes Modules
from ..core.config import global_conf
from ..core.job import Job

def mem(in1fastq, in2fastq=None, out_sam=None, read_group=None, ref=None, ini_section='bwa_kit'):
	other_options = global_conf.global_get(ini_section, 'other_options', required=False)
	out_bam = out_sam + ".aln.bam"
 
	return Job(
		[in1fastq, in2fastq, ref + ".bwt" if ref else None],
		[out_bam],
		[
			["bwa_kit", "module_bwakit"]
		],
        command="""\
bash <(run-bwamem {other_options}{read_group} \\
  {idxbase} \\
  {in1fastq}{in2fastq}{out_bam})""".format(
        other_options=" \\\n  " + other_options if other_options else "",
        read_group=" \\\n  -R " + read_group if read_group else "",
        idxbase=ref if ref else global_conf.global_get(ini_section, 'genome_bwa_index', param_type='filepath'),
        in1fastq=in1fastq,
        in2fastq=" \\\n  " + in2fastq if in2fastq else "",
        out_bam=" \\\n  -o " + out_bam if out_bam else ""
        ),
        removable_files=[out_bam]
    )

def bwa_postalt(input, output):
	return Job(
		[input],
		[output],
		[
			["bwa_kit", "module_bwakit"]
		],
		command="""\
$BWAKIT_PATH/k8 $BWAKIT_PATH/bwa-postalt.js {alt} \\
	{input} \\
	{output}""".format(
		alt=global_conf.global_get('bwa_kit', 'genome_bwa_alt_index', param_type='filepath'),
		input=input if input else "",
		output=output if output else ""
		),
		removable_files=[output]
	)
