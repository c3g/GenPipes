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
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def ensemble(lumpy, manta, cnvkit, wham, delly, gatk, bam, sample_name, workdir, outdir, isize_mean, isize_sd, output_vcf, breakseq=None):
    return Job(
        [lumpy, manta, cnvkit, wham, delly, breakseq, gatk, bam],
        [output_vcf],
        [
            ['metasv_ensemble', 'module_spades'],
            ['metasv_ensemble', 'module_age'],
            ['metasv_ensemble', 'module_python'],
            ['metasv_ensemble', 'module_bedtools'],
        ],
        command="""\
run_metasv.py {options} \\
  --num_threads {threads} \\
  --reference {genome} \\
  --lumpy_vcf {lumpy} \\
  --manta_vcf {manta} \\
  {cnvkit} \\
  {wham} \\
  {delly} \\
  {gatk} \\
  {breakseq} \\
  --bam {bam} \\
  --sample {sample_name} \\
  --spades spades.py \\
  --age age_align \\
  --workdir {workdir} \\
  --outdir {outdir} \\
  {isize_mean} \\
  {isize_sd}""".format(
        options=global_config_parser.param('metasv_ensemble', 'options'),
        threads=global_config_parser.param('metasv_ensemble', 'threads'),
        genome=global_config_parser.param('metasv_ensemble', 'genome_fasta', type='filepath'),
        lumpy=lumpy,
        manta=manta,
        wham="--wham_vcf " + wham if wham else "",
        delly="--delly_vcf " + delly if delly else "",
        gatk="--gatk_vcf " + gatk if gatk else "",
        cnvkit="--cnvkit_vcf " + cnvkit if cnvkit else "",
        breakseq="--breakseq_vcf " + breakseq if breakseq else "",
        bam=bam,
        sample_name=sample_name,
        workdir=workdir,
        outdir=outdir,
        isize_mean="--isize_mean " + isize_mean if isize_mean else "",
        isize_sd="--isize_sd " + isize_sd if isize_sd else "",
        )
    )
