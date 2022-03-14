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

def check(input_vcf, output_log, chr):
    output_dep = output_log + ".snp.strand.exclude"
    return Job(
        [input_vcf],
        [output_dep],
        [
            ['DEFAULT', 'module_shapeit'],
        ],
        command="""\
shapeit -check {options} \\
  --thread {threads} \\
  --input-vcf {input_vcf} \\
  -M {genome_map} \\
  --input-ref {hap_file} \\
  {legend_file} \\
  {sample_file} \\
  --output-log {output_log} ; true""".format(
        options=global_config_parser.param('shapeit', 'check_options'),
        threads=global_config_parser.param('shapeit', 'check_threads'),
        input_vcf=input_vcf,
        genome_map="${SHAPEIT_PATH}/ALL.integrated_phase1_SHAPEIT_16-06-14.nosing/genetic_map_chr" + chr + "_combined_b37.txt",
        hap_file="${SHAPEIT_PATH}/ALL.integrated_phase1_SHAPEIT_16-06-14.nosing/ALL.chr" + chr + ".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.haplotypes.gz",
        legend_file="${SHAPEIT_PATH}/ALL.integrated_phase1_SHAPEIT_16-06-14.nosing/ALL.chr" + chr + ".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz",
        sample_file="${SHAPEIT_PATH}/ALL.integrated_phase1_SHAPEIT_16-06-14.nosing/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample",
        output_log=output_log,
        )
    )

def phase(input_vcf, exclude_snps, output, output_log, chr):
    output_dep = output + ".haps"
    return Job(
        [input_vcf, exclude_snps],
        [output_dep],
        [
            ['DEFAULT', 'module_shapeit'],
        ],
        command="""\
shapeit {options} \\
  --thread {threads} \\
  --input-vcf {input_vcf} \\
  -M {genome_map} \\
  --input-ref {hap_file} \\
  {legend_file} \\
  {sample_file} \\
  --exclude-snp {exclude_snps} \\
  --output-log {output_log} \\
  -O {output}""".format(
        options=global_config_parser.param('shapeit', 'phase_options'),
        threads=global_config_parser.param('shapeit', 'phase_threads'),
        input_vcf=input_vcf,
        genome_map="${SHAPEIT_PATH}/ALL.integrated_phase1_SHAPEIT_16-06-14.nosing/genetic_map_chr" + chr + "_combined_b37.txt",
        hap_file="${SHAPEIT_PATH}/ALL.integrated_phase1_SHAPEIT_16-06-14.nosing/ALL.chr" + chr + ".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.haplotypes.gz",
        legend_file="${SHAPEIT_PATH}/ALL.integrated_phase1_SHAPEIT_16-06-14.nosing/ALL.chr" + chr + ".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz",
        sample_file="${SHAPEIT_PATH}/ALL.integrated_phase1_SHAPEIT_16-06-14.nosing/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample",
        exclude_snps=exclude_snps,
        output_log=output_log,
        output=output,
        )
    )
