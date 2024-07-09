################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def run(
    input_bam,
    sample,
    output_dir,
    junction_file=None,
    ini_section='sv_prep'
    ):
    inputs = [input_bam]
    if junction_file:
        inputs.append(junction_file)
    return Job(
        inputs,
        [
            os.path.join(output_dir, f"{sample}.sv_prep.bam"),
            os.path.join(output_dir, f"{sample}.sv_prep.junctions.csv")
        ],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_sv_prep']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $SVPREP_JAR  \\
  -threads {threads} \\
  -sample {sample} \\
  -bam_file {bam_file} \\
  -ref_genome {reference_sequence} \\
  -ref_genome_version {build} \\
  -blacklist_bed {blacklist_bed} \\
  -known_fusion_bed {known_fusion} \\
  {junction_file} \\
  -output_dir {outdir}""".format(
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
            ram=global_conf.global_get(ini_section, 'ram'),
            threads=global_conf.global_get(ini_section, 'threads'),
            sample=sample,
            bam_file=input_bam,
            reference_sequence=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
            build=global_conf.global_get(ini_section, 'assembly_alias2'),
            blacklist_bed=global_conf.global_get(ini_section, 'blacklist_bed', param_type='filepath'),
            known_fusion=global_conf.global_get(ini_section, 'known_fusion', param_type='filepath'),
            junction_file=f"-existing_junction_file {junction_file}" if junction_file else "",
            outdir=output_dir
        )
    )
