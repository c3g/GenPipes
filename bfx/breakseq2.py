################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
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

from core.config import *
from core.job import *

def run(input, sample_name, output_dir):
	output = os.path.join(output_dir, "breakseq.vcf.gz")
	return Job(
        [input],
        [output],
        [
            ['run_breakseq2', 'module_python'],
	        ['run_breakseq2', 'module_breakseq2'],
	        ['run_breakseq2', 'module_samtools'],
	        ['run_breakseq2', 'module_bwa']
        ],
        command="""\
run_breakseq2.py {options} --bwa bwa --samtools samtools \\
    --nthreads {threads} \\
    --reference {genome}  \\
    --bplib_gff {gff} \\
    --bams {input} \\
    --sample {sample} \\
    --work {output}""".format(
            options=config.param('run_breakseq2','options'),
            threads=config.param('run_breakseq2','threads'),
            genome=config.param('run_breakseq2','genome_fasta',type='filepath'),
	        gff=config.param('run_breakseq2','gff',type='filepath'),
	        output=output_dir,
	        input=input,
	        sample=sample_name,
        )
    )
