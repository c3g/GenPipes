################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job


def run(
    input_bam,
    output,
    tmp_dir,
    region=None,
    ini_section='deepvariant'
    ):
    """
    Variant calling with DeepVariant.

    :return: a job for DeepVariant variant calling.
    """

    return Job(
        [input_bam],
        [output],
        [
            [ini_section, "module_deepvariant"],
            [ini_section, "module_apptainer"]
        ],
        command="""\
apptainer exec \\
  --bind /cvmfs/soft.mugqic,/cvmfs/ref.mugqic,{extra_paths_to_bind} \\
  -e $DEEPVARIANT_ROOT \\
  run_deepvariant {other_options} \\
  --reads {input_bam} \\
  --ref {genome_fasta} \\
  --model_type {model_type} \\
  {region} \\
  --num_shards {threads} \\
  --intermediate_results_dir {tmp_dir} \\
  --output_vcf {output}""".format(
            extra_paths_to_bind=global_conf.global_get(ini_section, 'extra_paths_to_bind', required=False),
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            input_bam=input_bam,
            genome_fasta=global_conf.global_get(ini_section, 'genome_fasta'),
            model_type=global_conf.global_get(ini_section, 'model_type'),
            threads=global_conf.global_get(ini_section, 'threads'),
            tmp_dir=tmp_dir,
            region="--regions " + region if region else "",
            output=output
        ),
        removable_files=tmp_dir
    )
