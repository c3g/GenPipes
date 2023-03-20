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
import logging
import os

# MUGQIC Modules
from core.config import config
from core.job import Job

def filter(
    gridss_vcf,
    gripss_vcf,
    gripss_filtered_vcf,
    output_id,
    sample,
    reference
    ):
    return Job(
        [gridss_vcf],
        [
            gripss_vcf,
            gripss_filtered_vcf
        ],
        [
            ['gripss_filter', 'module_java'],
            ['gripss_filter', 'module_gripss']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $GRIPSS_JAR  \\
  -ref_genome {reference_sequence} \\
  -known_hotspot_file {known_hotspot} \\
  -pon_sgl_file {pon_sgl} \\
  -pon_sv_file {pon_sv} \\
  -vcf {gridss_vcf} \\
  -output_id {output_id} \\
  -output_dir {outdir} \\
  -sample {sample} \\
  -referecence {reference}""".format(
            tmp_dir=config.param('gripss_filter', 'tmp_dir'),
            java_other_options=config.param('gripss_filter', 'java_other_options'),
            ram=config.param('gripss_filter', 'ram'),
            reference_sequence=config.param('gripss_filter', 'genome_fasta', param_type='filepath'),
            known_hotspot=config.param('gripss_filter', 'known_hotspot', param_type='filepath'),
            pon_sgl=config.param('gripss_filter', 'pon_sgl', param_type='filepath'),
            pon_sv=config.param('gripss_filter', 'pon_sv', param_type='filepath'),
            gridss_vcf=gridss_vcf,
            output_id=output_id,
            outdir=os.path.dirname(gripss_vcf),
            sample=sample,
            reference=reference
        )
    )
