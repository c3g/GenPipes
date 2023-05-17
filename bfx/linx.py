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

def somatic(
    purple_vcf,
    purple_dir,
    tumor_name,
    output_dir
    ):
    purple_input = os.path.join(purple_dir, tumor_name + ".purple.purity.tsv")
    linx_output = os.path.join(output_dir, tumor_name + ".linx.vis_sv_data.tsv")    

    return Job(
        [purple_vcf, purple_input],
        [linx_output],
        [
            ['linx', 'module_java'],
            ['linx', 'module_R'],
            ['linx', 'module_perl'],
            ['linx', 'module_circos'],
            ['linx', 'module_linx']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $LINX_JAR \\
  -threads {threads} \\
  -ref_genome_version {build} \\
  -sample {sample} \\
  -purple_dir {purple_dir} \\
  -sv_vcf {purple_vcf} \\
  -output_dir {outdir} \\
  -fragile_site_file {fragile_site} \\
  -line_element_file {line_element} \\
  -ensembl_data_dir {ensembl_data} \\
  -check_fusions \\
  -known_fusion_file {known_fusion} \\
  -driver_gene_panel {driver_gene} \\
  -check_drivers \\
  -write_vis_data """.format(
            tmp_dir=config.param('linx', 'tmp_dir'),
            java_other_options=config.param('linx', 'java_other_options'),
            ram=config.param('linx', 'ram'),
            threads=config.param('linx', 'threads'),
            build=config.param('linx', 'assembly_alias2'),
            sample=tumor_name,
            purple_dir=purple_dir,
            purple_vcf=purple_vcf,
            outdir=output_dir,
            fragile_site=config.param('linx', 'fragile_site', param_type='filepath'),
            line_element=config.param('linx', 'line_element', param_type='filepath'),
            ensembl_data=config.param('linx', 'ensembl_data', param_type='dirpath'),
            known_fusion=config.param('linx', 'known_fusion', param_type='filepath'),
            driver_gene=config.param('linx', 'driver_gene', param_type='filepath')
        )
    )

def germline(
    purple_vcf,
    tumor_name,
    output_dir
    ):
    purple_input = os.path.join(os.path.dirname(purple_vcf), tumor_name + ".purple.purity.tsv")
    linx_output = os.path.join(output_dir, tumor_name + ".linx.vis_sv_data.tsv")

    return Job(
        [purple_vcf, purple_input],
        [linx_output],
        [
            ['linx', 'module_java'],
            ['linx', 'module_R'],
            ['linx', 'module_perl'],
            ['linx', 'module_circos'],
            ['linx', 'module_linx']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $LINX_JAR \\
  -threads {threads} \\
  -ref_genome_version {build} \\
  -sample {sample} \\
  -germline \\
  -sv_vcf {purple_vcf} \\
  -output_dir {outdir} \\
  -fragile_site_file {fragile_site} \\
  -line_element_file {line_element} \\
  -ensembl_data_dir {ensembl_data} \\
  -driver_gene_panel {driver_gene}""".format(
            tmp_dir=config.param('linx', 'tmp_dir'),
            java_other_options=config.param('linx', 'java_other_options'),
            ram=config.param('linx', 'ram'),
            threads=config.param('linx', 'threads'),
            build=config.param('linx', 'assembly_alias2'),
            sample=tumor_name,
            purple_vcf=purple_vcf,
            outdir=output_dir,
            fragile_site=config.param('linx', 'fragile_site', param_type='filepath'),
            line_element=config.param('linx', 'line_element', param_type='filepath'),
            ensembl_data=config.param('linx', 'ensembl_data', param_type='dirpath'),
            driver_gene=config.param('linx', 'driver_gene', param_type='filepath')
        )
    )


def plot(
    sample_name,
    linx_dir
    ):
    linx_input = os.path.join(linx_dir, sample_name + ".linx.vis_sv_data.tsv")

    return Job(
        [linx_input],
        [
            os.path.join(linx_dir, "plot"),
            os.path.join(linx_dir, "circos")
        ],
        [
            ['linx', 'module_java'],
            ['linx', 'module_R'],
            ['linx', 'module_perl'],
            ['linx', 'module_circos'],
            ['linx', 'module_linx']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -cp $LINX_JAR com.hartwig.hmftools.linx.visualiser.SvVisualiser \\
  -threads {threads} \\
  -sample {sample} \\
  -plot_out {linx_dir}/plot/ \\
  -data_out {linx_dir}/circos/ \\
  -vis_file_dir {linx_dir} \\
  -circos circos""".format(
            tmp_dir=config.param('linx', 'tmp_dir'),
            java_other_options=config.param('linx', 'java_other_options'),
            ram=config.param('linx', 'ram'),
            threads=config.param('linx', 'threads'),
            sample=sample_name,
            linx_dir=linx_dir
        )
    )
