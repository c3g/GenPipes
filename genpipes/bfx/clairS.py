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

import os

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job


def run(
    normal_bam,
    tumor_bam,
    output_dir,
    tumor_name,
    platform,
    region=None,
    ini_section='clairS'
    ):
    """
    Call somatic small variants with ClairS

    :return: a job for ClairS somatic small variant calling
    """

    outputs = [
        os.path.join(output_dir, "output.vcf.gz"),
        os.path.join(output_dir, "clair3_normal_germline_output.vcf.gz"),
        os.path.join(output_dir, "indel.vcf.gz"),
        os.path.join(output_dir, "snv.vcf.gz")
    ]

    return Job(
        [normal_bam, tumor_bam],
        outputs,
        [
            [ini_section, "module_clairS"],
        ],
        command="""\
run_clairs {other_options} \\
  --normal_bam_fn={normal_bam} \\
  --tumor_bam_fn {tumor_bam} \\
  --ref_fn={genome_fasta} \\
  --threads={threads} \\
  --platform="{platform}" \\
  --clair3_model_path={model_path} \\
  {region} {sites_to_call} \\
  --sample_name={tumor_name} \\
  --output_dir={output}""".format(
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            normal_bam=normal_bam,
            tumor_bam=tumor_bam,
            genome_fasta=global_conf.global_get(ini_section, 'genome_fasta'),
            threads=global_conf.global_get(ini_section, 'threads'),
            platform=platform,
            model_path=global_conf.global_get(ini_section, 'model_path', param_type='dirpath'),
            region=region if region else "",
            sites_to_call="--vcf_fn=" + global_conf.global_get(ini_section, 'sites_to_call', required=False, param_type='filepath') if global_conf.global_get(ini_section, 'sites_to_call', required=False) else "",
            tumor_name=tumor_name,
            output=output_dir
        )
    )
