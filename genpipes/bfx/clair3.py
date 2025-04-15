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
    input_bam,
    output_dir,
    sample_name,
    platform,
    region=None,
    ini_section='clair3'
    ):
    """
    Call variants with Clair3

    :return: a job for Clair3 variant calling
    """

    outputs = [
        os.path.join(output_dir, "pileup.vcf.gz"),
        os.path.join(output_dir, "full_alignment.vcf.gz"),
        os.path.join(output_dir, "merge_output.vcf.gz"),
        os.path.join(output_dir, "phased_merge_output.vcf.gz")
    ]

    return Job(
        [input_bam],
        outputs,
        [
            [ini_section, "module_clair3"],
        ],
        command="""\
run_clair3.sh {other_options} \\
  --bam_fn={input_bam} \\
  --ref_fn={genome_fasta} \\
  --threads={threads} \\
  --platform="{platform}" \\
  --model_path={model_path} \\
  {region} {sites_to_call} \\
  --sample_name={sample_name} \\
  --output={output}""".format(
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            input_bam=input_bam,
            genome_fasta=global_conf.global_get(ini_section, 'genome_fasta'),
            threads=global_conf.global_get(ini_section, 'threads'),
            platform=platform,
            model_path=global_conf.global_get(ini_section, 'model_path', param_type='dirpath'),
            region=region if region else "",
            sites_to_call="--vcf_fn=" + global_conf.global_get(ini_section, 'sites_to_call', required=False, param_type='filepath') if global_conf.global_get(ini_section, 'sites_to_call', required=False) else "",
            sample_name=sample_name,
            output=output_dir
        )
    )
