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
import logging
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def run(
        somatic_snv,
        tumor_name,
        output_dir,
        output_file,
        ini_section='pave_annotate'
    ):

    return Job(
        [somatic_snv],
        [output_file],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_R'],
            [ini_section, 'module_pave'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PAVE_JAR \\
  {other_options} \\
  -threads {threads} \\
  -ref_genome_version {reference_sequence_version} \\
  -ref_genome {reference_sequence} \\
  -sample {tumor_name} \\
  -vcf_file {somatic_snv} \\
  -ensembl_data_dir {ensembl_data_dir} \\
  -driver_gene_panel {driver_gene_panel} \\
  -gnomad_freq_dir {gnomad_freq_dir} \\
  -mappability_bed {mappability_bed} \\
  -clinvar_vcf {clinvar_vcf} \\
  -output_dir {output_dir} \\
  -output_vcf_file {output_file}""".format(
            tmp_dir=config.param(ini_section, 'tmp_dir'),
            java_other_options=config.param(ini_section, 'java_other_options'),
            other_options=config.param(ini_section, 'other_options', required=False),
            threads=config.param(ini_section, 'threads'),
            ram=config.param(ini_section, 'ram'),
            reference_sequence_version=config.param('DEFAULT', 'assembly_alias2'),
            reference_sequence=config.param(ini_section, 'genome_fasta', param_type='filepath'),
            tumor_name=tumor_name,
            somatic_snv=somatic_snv,
            ensembl_data_dir=config.param(ini_section, 'ensembl_data_dir', param_type='dirpath'),
            driver_gene_panel=config.param(ini_section, 'driver_gene_panel', param_type='filepath'),
            gnomad_freq_dir=config.param(ini_section, 'gnomad_freq_dir', param_type='dirpath'),
            mappability_bed=config.param(ini_section, 'mappability_bed', param_type='filepath'),
            clinvar_vcf=config.param(ini_section, 'clinvar_vcf', param_type='filepath'),
            output_dir=output_dir,
            output_file=output_file
        )
    )