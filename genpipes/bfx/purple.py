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
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import os

# GenPipes Modules
from ..core.config import global_conf
from ..core.job import Job

def run(
    amber,
    cobalt,
    normal_name,
    tumor_name,
    output_dir,
    ensembl_data_dir,
    somatic_snv,
    structural_sv=None,
    sv_recovery=None,
    somatic_hotspots=None,
    germline_hotspots=None,
    driver_gene_panel=None,
    ini_section='purple'
    ):

    amber_input = os.path.join(amber, tumor_name + ".amber.baf.pcf")
    cobalt_input = os.path.join(cobalt, tumor_name + ".cobalt.ratio.pcf")
    input_files = [amber_input, cobalt_input, somatic_snv]

    purple_outputs = [
        os.path.join(output_dir, tumor_name + ".purple.purity.tsv"),
        os.path.join(output_dir, tumor_name + ".purple.qc"),
        os.path.join(output_dir, tumor_name + ".purple.driver.catalog.somatic.tsv"),
        os.path.join(output_dir, tumor_name + ".purple.driver.catalog.germline.tsv")
    ]

    if structural_sv is not None and sv_recovery is not None:
        input_files.append(structural_sv)
        input_files.append(sv_recovery)
        purple_sv = os.path.join(output_dir, tumor_name + ".purple.sv.vcf.gz")
        driver_somatic = os.path.join(output_dir, tumor_name + ".purple.driver.catalog.somatic.tsv")
        driver_germline = os.path.join(output_dir, tumor_name + ".purple.driver.catalog.germline.tsv")
        purple_outputs.append(purple_sv)
        purple_outputs.append(driver_somatic)
        purple_outputs.append(driver_germline)

    if structural_sv is not None:
        circos_plot = os.path.join(output_dir, "plot", tumor_name + ".circos.png")
        purple_outputs.append(circos_plot)

    return Job(
        input_files,
        purple_outputs,
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_R'],
            [ini_section, 'module_perl'],
            [ini_section, 'module_circos'],
            [ini_section, 'module_cobalt'],
            [ini_section, 'module_purple'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PURPLE_JAR \\
  {other_options} \\
  -threads {threads} \\
  -ref_genome_version {reference_sequence_version} \\
  -ref_genome {reference_sequence} \\
  -reference {reference} \\
  -tumor {tumor} \\
  -cobalt {cobalt} \\
  -amber {amber} \\
  -gc_profile {gc_profile} \\
  -ensembl_data_dir {ensembl_data_dir} \\
  -somatic_vcf {somatic_snv} {structural_sv} {sv_recovery} {somatic_hotspots} {germline_hotspots} {driver_gene_panel} {circos} \\
  -output_dir {output_dir}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            run_drivers="\\\n  -run_drivers" if structural_sv else "",
            threads=global_conf.global_get(ini_section, 'threads'),
            ram=global_conf.global_get(ini_section, 'ram'),
            gc_profile=global_conf.global_get(ini_section, 'gc_profile'),
            reference_sequence_version=global_conf.global_get('DEFAULT', 'assembly_alias2'),
            reference_sequence=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
            reference=normal_name,
            tumor=tumor_name,
            amber=amber,
            cobalt=cobalt,
            ensembl_data_dir=ensembl_data_dir,
            somatic_snv=somatic_snv,
            structural_sv=" \\\n  -somatic_sv_vcf " +  structural_sv if structural_sv else "",
            sv_recovery=" \\\n  -sv_recovery_vcf " +  sv_recovery if sv_recovery else "",
            somatic_hotspots=" \\\n  -somatic_hotspots " + somatic_hotspots if somatic_hotspots else "",
            germline_hotspots=" \\\n  -germline_hotspots " + germline_hotspots if germline_hotspots else "",
            driver_gene_panel=" \\\n  -driver_gene_panel " + driver_gene_panel if driver_gene_panel else "",
            circos="\\\n  -circos `which circos`" if structural_sv else "",   # full path to circos binary required
            output_dir=output_dir
        )
    )

def strelka2_convert(
        input,
        output,
        ini_section='purple_convert_strelka2'
):

    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_R'],
            [ini_section, 'module_purple'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -cp $PURPLE_JAR com.hartwig.hmftools.purple.tools.AnnotateStrelkaWithAllelicDepth \\
  -in {input} \\
  -out {output}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            input=input,
            output=output,
        )
    )

def parse_purity_metrics_pt(input_file):
    """
    Parse purity metrics from Purple output.
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export purity=`awk 'NR>1{{printf "%.0f", $1*100}}' {input_file}`"""
        )
