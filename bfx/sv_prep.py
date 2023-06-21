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
from core.config import config
from core.job import Job

def run(
    input_bam,
    sample,
    output_dir,
    junction_file=None
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
            ['sv_prep', 'module_java'],
            ['sv_prep', 'module_sv_prep']
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
            tmp_dir=config.param('sv_prep', 'tmp_dir'),
            java_other_options=config.param('sv_prep', 'java_other_options'),
            ram=config.param('sv_prep', 'ram'),
            threads=config.param('sv_prep', 'threads'),
            sample=sample,
            bam_file=input_bam,
            reference_sequence=config.param('sv_prep', 'genome_fasta', param_type='filepath'),
            build=config.param('sv_prep', 'assembly_alias2'),
            blacklist_bed=config.param('sv_prep', 'blacklist_bed', param_type='filepath'),
            known_fusion=config.param('sv_prep', 'known_fusion', param_type='filepath'),
            junction_file=f"-existing_junction_file {junction_file}" if junction_file else "",
            outdir=output_dir
        )
    )
