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
    output_dir,
    ini_section='linx_somatic'
    ):
    purple_input = os.path.join(purple_dir, tumor_name + ".purple.purity.tsv")
    linx_output = [
            os.path.join(output_dir, tumor_name + ".linx.breakend.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.clusters.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.driver.catalog.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.drivers.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.fusion.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.links.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.svs.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.vis_copy_number.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.vis_fusion.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.vis_gene_exon.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.vis_protein_domain.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.vis_segments.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.vis_sv_data.tsv")
        ]

    return Job(
        [purple_vcf, purple_input],
        linx_output,
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_R'],
            [ini_section, 'module_perl'],
            [ini_section, 'module_circos'],
            [ini_section, 'module_linx']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $LINX_JAR \\
  {options} \\
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
            tmp_dir=config.param(ini_section, 'tmp_dir'),
            java_other_options=config.param(ini_section, 'java_other_options'),
            ram=config.param(ini_section, 'ram'),
            options=config.param(ini_section, 'options'),
            threads=config.param(ini_section, 'threads'),
            build=config.param(ini_section, 'assembly_alias2'),
            sample=tumor_name,
            purple_dir=purple_dir,
            purple_vcf=purple_vcf,
            outdir=output_dir,
            fragile_site=config.param(ini_section, 'fragile_site', param_type='filepath'),
            line_element=config.param(ini_section, 'line_element', param_type='filepath'),
            ensembl_data=config.param(ini_section, 'ensembl_data', param_type='dirpath'),
            known_fusion=config.param(ini_section, 'known_fusion', param_type='filepath'),
            driver_gene=config.param(ini_section, 'driver_gene', param_type='filepath')
        )
    )

def germline(
    purple_vcf,
    tumor_name,
    output_dir,
    ini_section='linx_germline'
    ):
    purple_input = os.path.join(os.path.dirname(purple_vcf), tumor_name + ".purple.purity.tsv")
    linx_output = [
            os.path.join(output_dir, tumor_name + ".linx.germline.clusters.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.germline.disruption.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.germline.driver.catalog.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.germline.links.tsv"),
            os.path.join(output_dir, tumor_name + ".linx.germline.svs.tsv")
        ]

    return Job(
        [purple_vcf, purple_input],
        linx_output,
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_R'],
            [ini_section, 'module_perl'],
            [ini_section, 'module_circos'],
            [ini_section, 'module_linx']
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
            tmp_dir=config.param(ini_section, 'tmp_dir'),
            java_other_options=config.param(ini_section, 'java_other_options'),
            ram=config.param(ini_section, 'ram'),
            threads=config.param(ini_section, 'threads'),
            build=config.param(ini_section, 'assembly_alias2'),
            sample=tumor_name,
            purple_vcf=purple_vcf,
            outdir=output_dir,
            fragile_site=config.param(ini_section, 'fragile_site', param_type='filepath'),
            line_element=config.param(ini_section, 'line_element', param_type='filepath'),
            ensembl_data=config.param(ini_section, 'ensembl_data', param_type='dirpath'),
            driver_gene=config.param(ini_section, 'driver_gene', param_type='filepath')
        )
    )


def plot(
    sample_name,
    linx_dir,
    ini_section='linx_plot'
    ):
    linx_input = os.path.join(linx_dir, sample_name + ".linx.vis_sv_data.tsv")

    return Job(
        [linx_input],
        [
            os.path.join(linx_dir, "plot"),
            os.path.join(linx_dir, "circos")
        ],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_R'],
            [ini_section, 'module_perl'],
            [ini_section, 'module_circos'],
            [ini_section, 'module_linx']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -cp $LINX_JAR com.hartwig.hmftools.linx.visualiser.SvVisualiser \\
  -threads {threads} \\
  -sample {sample} \\
  -plot_out {linx_dir}/plot/ \\
  -data_out {linx_dir}/circos/ \\
  -vis_file_dir {linx_dir} \\
  -circos circos""".format(
            tmp_dir=config.param(ini_section, 'tmp_dir'),
            java_other_options=config.param(ini_section, 'java_other_options'),
            ram=config.param(ini_section, 'ram'),
            threads=config.param(ini_section, 'threads'),
            sample=sample_name,
            linx_dir=linx_dir
        )
    )
