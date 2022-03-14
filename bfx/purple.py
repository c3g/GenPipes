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

def run(amber, cobalt, normal_name, tumor_name, output_dir, somatic_snv):
    amber_input = os.path.join(amber, tumor_name + ".amber.baf.pcf")
    cobalt_input = os.path.join(cobalt, tumor_name + ".cobalt.ratio.pcf")
    purple_output = os.path.join(output_dir, tumor_name + ".purple.purity.tsv")
    
    return Job(
        [amber_input, cobalt_input, somatic_snv],
        [purple_output],
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
  -ref_genome {reference_sequence} \\
  -reference {reference} \\
  -tumor {tumor} \\
  -cobalt {cobalt} \\
  -amber {amber} \\
  -gc_profile {gc_profile} \\
  -circos circos \\
  -somatic_vcf {somatic_snv} \\
  -output_dir {output_dir}""".format(
        tmp_dir=global_config_parser.param('purple', 'tmp_dir'),
        java_other_options=global_config_parser.param('purple', 'java_other_options'),
        ram=global_config_parser.param('purple', 'ram'),
        threads=global_config_parser.param('purple', 'threads'),
        gc_profile=global_config_parser.param('purple', 'gc_profile'),
        reference_sequence=global_config_parser.param('purple', 'genome_fasta', param_type='filepath'),
        reference=normal_name,
        tumor=tumor_name,
        amber=amber,
        cobalt=cobalt,
        somatic_snv=somatic_snv,
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
        tmp_dir=global_config_parser.param('purple', 'tmp_dir'),
        java_other_options=global_config_parser.param('purple', 'java_other_options'),
        ram=global_config_parser.param('purple', 'ram'),
        input=input,
        output=output,
        )
    )
