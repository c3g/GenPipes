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
from core.config import *
from core.job import *

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
    driver_gene_panel=None
    ):
    amber_input = os.path.join(amber, tumor_name + ".amber.baf.pcf")
    cobalt_input = os.path.join(cobalt, tumor_name + ".cobalt.ratio.pcf")
<<<<<<< HEAD
    purple_outputs = [
        os.path.join(output_dir, tumor_name + ".purple.purity.tsv"),
        os.path.join(output_dir, tumor_name + ".purple.qc")
    ]
||||||| parent of c2eac9fa (Tumor Pair - SV : revamping SV protocol of Tumor Pair with use of amber, purple, gridss, etc...)
    purple_output = os.path.join(output_dir, tumor_name + ".purple.purity.tsv")
=======
    input_files = [amber_input, cobalt_input, somatic_snv]

    purple_output = [os.path.join(output_dir, tumor_name + ".purple.purity.tsv")]
>>>>>>> c2eac9fa (Tumor Pair - SV : revamping SV protocol of Tumor Pair with use of amber, purple, gridss, etc...)
    

    if structural_sv is not None and sv_recovery is not None:
        input_files.append(structural_sv)
        input_files.append(sv_recovery)
        purple_sv = os.path.join(output_dir, tumor_name + ".purple.sv.vcf.gz")
        purple_output.append(purple_sv)    

    return Job(
<<<<<<< HEAD
        [amber_input, cobalt_input, somatic_snv],
        purple_outputs,
||||||| parent of c2eac9fa (Tumor Pair - SV : revamping SV protocol of Tumor Pair with use of amber, purple, gridss, etc...)
        [amber_input, cobalt_input, somatic_snv],
        [purple_output],
=======
        input_files,
        purple_output,
>>>>>>> c2eac9fa (Tumor Pair - SV : revamping SV protocol of Tumor Pair with use of amber, purple, gridss, etc...)
        [
            ['purple', 'module_java'],
            ['purple', 'module_R'],
            ['purple', 'module_perl'],
            ['purple', 'module_circos'],
            ['purple', 'module_cobalt'],
            ['purple', 'module_purple'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PURPLE_JAR \\
  {other_options} \\
  {run_drivers} \\
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
            tmp_dir=config.param('purple', 'tmp_dir'),
            java_other_options=config.param('purple', 'java_other_options'),
            other_options=config.param('purple', 'other_options', required=False),
            run_drivers="\\\n  -run_drivers" if structural_sv else "",
            threads=config.param('purple', 'threads'),
            ram=config.param('purple', 'ram'),
            gc_profile=config.param('purple', 'gc_profile'),
            reference_sequence_version=config.param('DEFAULT', 'assembly_alias2'),
            reference_sequence=config.param('purple', 'genome_fasta', param_type='filepath'),
            reference=normal_name,
            tumor=tumor_name,
            amber=amber,
            cobalt=cobalt,
            ensembl_data_dir=ensembl_data_dir,
            somatic_snv=somatic_snv,
            structural_sv=" \\\n  -structural_vcf " +  structural_sv if structural_sv else "",
            sv_recovery=" \\\n  -sv_recovery_vcf " +  sv_recovery if sv_recovery else "",
            somatic_hotspots=" \\\n  -somatic_hotspots " + somatic_hotspots if somatic_hotspots else "",
            germline_hotspots=" \\\n  -germline_hotspots " + germline_hotspots if germline_hotspots else "",
            driver_gene_panel=" \\\n  -driver_gene_panel " + driver_gene_panel if driver_gene_panel else "",
            circos="\\\n  -circos circos" if structural_sv else "", 
            output_dir=output_dir
        )
    )

def strelka2_convert(input, output):

    return Job(
        [input],
        [output],
        [
            ['purple', 'module_java'],
            ['purple', 'module_R'],
            ['purple', 'module_purple'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -cp $PURPLE_JAR com.hartwig.hmftools.purple.tools.AnnotateStrelkaWithAllelicDepth \\
  -in {input} \\
  -out {output}""".format(
            tmp_dir=config.param('purple', 'tmp_dir'),
            java_other_options=config.param('purple', 'java_other_options'),
            ram=config.param('purple', 'ram'),
            input=input,
            output=output,
        )
    )
