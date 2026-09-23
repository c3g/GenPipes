################################################################################
# Copyright (C) 2026 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
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
    tumor_name,
    snv_indel_vcf,
    sv_vcf,
    output_dir,
    ini_section='chord'
    ):

    chord_outputs = [
        os.path.join(output_dir, tumor_name + ".chord.mutation_contexts.tsv"),
        os.path.join(output_dir, tumor_name + ".chord.prediction.tsv")
    ]

    return Job(
        [snv_indel_vcf, sv_vcf],
        chord_outputs,
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_R'],
            [ini_section, 'module_chord']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $CHORD_JAR \\
  {other_options} \\
  -threads {threads} \\
  -ref_genome {reference_sequence} \\
  -sample {tumor_name} \\
  -snv_indel_vcf_file {snv_indel_vcf} \\
  -sv_vcf_file {sv_vcf} \\
  -output_dir {output_dir}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            threads=global_conf.global_get(ini_section, 'threads'),
            ram=global_conf.global_get(ini_section, 'ram'),
            reference_sequence=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
            tumor_name=tumor_name,
            snv_indel_vcf=snv_indel_vcf,
            sv_vcf=sv_vcf,
            output_dir=output_dir
        )
    )
