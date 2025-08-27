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
    sample_name,
    germline_vcf = None,
    ini_section='savana'
    ):
    """
    Call somatic structural variants and CNAs with Savana

    :return: a job for Savana somatic SV and CNA calling
    """

    outputs = [
        os.path.join(output_dir, f"{sample_name}_segmented_absolute_copy_number.tsv"),
        os.path.join(output_dir, f"{sample_name}.classified.somatic.vcf"),
        os.path.join(output_dir, f"{sample_name}.classified.somatic.bedpe"),
        os.path.join(output_dir, f"{sample_name}.classified.vcf"),
        os.path.join(output_dir, f"{sample_name}_fitted_purity_ploidy.tsv"),
        os.path.join(output_dir, f"{sample_name}.inserted_sequences.fa"),
        os.path.join(output_dir, f"{sample_name}_raw_read_counts.tsv"),
        os.path.join(output_dir, f"{sample_name}.sv_breakpoints.vcf"),
        os.path.join(output_dir, f"{sample_name}.sv_breakpoints.bedpe"),
        os.path.join(output_dir, f"{sample_name}.sv_breakpoints_read_support.tsv"),
        os.path.join(output_dir, f"{sample_name}_ranked_solutions.tsv"),
        os.path.join(output_dir, f"{sample_name}_raw_read_counts.tsv"),
        os.path.join(output_dir, f"{sample_name}_read_counts_mnorm_log2r_segmented.tsv"),
        os.path.join(output_dir, f"{sample_name}_allele_counts_hetSNPs.bed")
    ]

    return Job(
        [normal_bam, tumor_bam, germline_vcf],
        outputs,
        [
            [ini_section, "module_savana"],
        ],
        command="""\
savana {other_options} \\
  --tumour {tumor_bam} \\
  --normal {normal_bam} \\
  --sample {sample_name} \\
  --ref {genome_fasta} \\
  --threads {threads} \\
  --cna_threads {cna_threads} \\
  {germline_vcf} \\
  {mapq} {cn_step_change} \\
  {contigs} \\
  --outdir {output}""".format(
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            normal_bam=normal_bam,
            tumor_bam=tumor_bam,
            sample_name=sample_name,
            genome_fasta=global_conf.global_get(ini_section, 'genome_fasta'),
            threads=global_conf.global_get(ini_section, 'threads'),
            cna_threads=global_conf.global_get(ini_section, 'cna_threads'),
            germline_vcf="--snp_vcf " + germline_vcf if germline_vcf else "",
            mapq="--mapq " + global_conf.global_get(ini_section, 'mapq'),
            cn_step_change=global_conf.global_get(ini_section, 'cn_step_change', required=False),
            contigs="--contigs " + global_conf.global_get(ini_section, 'contigs') if global_conf.global_get(ini_section, 'contigs') else "",            
            output=output_dir
        )
    )
