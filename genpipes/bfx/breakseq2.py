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

import os

from ..core.config import global_conf
from ..core.job import Job

def run(input, sample_name, output_dir, ini_section='run_breakseq2'):
	output = os.path.join(output_dir, "breakseq.vcf.gz")
	return Job(
        [input],
        [output],
        [
            [ini_section, 'module_python'],
	        [ini_section, 'module_breakseq2'],
	        [ini_section, 'module_samtools'],
	        [ini_section, 'module_bwa']
        ],
        command="""\
bwa_path=`which bwa`; \\
samtools_path=`which samtools`; \\
run_breakseq2.py {options} --bwa "$bwa_path" --samtools "$samtools_path" \\
    --nthreads {threads} \\
    --reference {genome}  \\
    --bplib_gff {gff} \\
    --bams {input} \\
    --sample {sample} \\
    --work {output}""".format(
            options=global_conf.global_get(ini_section,'options'),
            threads=global_conf.global_get(ini_section,'threads'),
            genome=global_conf.global_get(ini_section,'genome_fasta', param_type='filepath'),
	        gff=global_conf.global_get(ini_section,'gff', param_type='filepath'),
	        output=output_dir,
	        input=input,
	        sample=sample_name
        )
    )
