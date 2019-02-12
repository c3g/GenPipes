#!/usr/bin/env python

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

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def ensemble(input_lumpy, input_manta, input_cnvkit, input_wham, input_delly, input_gatk, input_bam, sample_name, workdir, outdir, isize_mean, isize_sd, output_vcf):
    return Job(
        [input_lumpy, input_manta, input_cnvkit, input_wham, input_delly, input_gatk, input_bam],
        [output_vcf],
        [
            #['metasv_ensemble', 'module_spades'],
            #['metasv_ensemble', 'module_age'],
            ['metasv_ensemble', 'module_python'],
            ['metasv_ensemble', 'module_bedtools'],
        ],
        command="""\
run_metasv.py {options} \\
  --num_threads {threads} \\
  --reference {genome} \\
  --lumpy_vcf {input_lumpy} \\
  --manta_vcf {input_manta} \\
  {input_cnvkit} \\
  {input_wham} \\
  {input_delly} \\
  {input_gatk} \\
  --bam {input_bam} \\
  --sample {sample_name} \\
  --spades $SPAdes_BIN/spades.py \\
  --age $AGE_BIN/age_align \\
  --workdir {workdir} \\
  --outdir {outdir} \\
  {isize_mean} \\
  {isize_sd}""".format(
        options=config.param('metasv_ensemble','options'),
        threads=config.param('metasv_ensemble','threads'),
        genome=config.param('metasv_ensemble','genome_fasta',type='filepath'),        
        input_lumpy=input_lumpy,
        input_manta=input_manta,
        input_wham="--wham_vcf " + input_wham if input_wham else "",
        input_delly="--delly_vcf " + input_delly if input_delly else "",
        input_gatk="--gatk_vcf " + input_gatk if input_gatk else "",
        input_cnvkit="--cnvkit_vcf " + input_cnvkit if input_cnvkit else "",
        input_bam=input_bam,
        sample_name=sample_name,
        workdir=workdir,
        outdir=outdir,
        isize_mean="--isize_mean " + isize_mean if isize_mean else "",
        isize_sd="--isize_sd " + isize_sd if isize_sd else "",
        )
    )
