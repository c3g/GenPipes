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


def annotate(
    input_vcf,
    output,
    input_sv_indel=None,
    ini_section='annotSV'
    ):
    """
    Annotate and rate structural variants with AnnotSV.

    :return: a job for AnnotSV annotation
    """

    return Job(
        [input_vcf, input_sv_indel],
        [output],
        [
            [ini_section, "module_annotsv"],
            [ini_section, "module_bedtools"],
            [ini_section, "module_bcftools"]

        ],
        command="""\
AnnotSV {other_options} \\
  -SVinputFile {input_vcf} \\
  {input_sv_indel} \\
  -annotationsDir $ANNOTSV/share/AnnotSV \\
  -genomeBuild {genome_build} \\
  {candidate_genes} \\
  -outputFile {output}""".format(
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            input_vcf=input_vcf,
            input_sv_indel="-snvIndelFiles " + input_sv_indel if input_sv_indel else "",
            genome_build=global_conf.global_get(ini_section, 'assembly'),
            candidate_genes="-candidateGenesFile " + global_conf.global_get(ini_section, 'candidate_genes', required=False) if global_conf.global_get(ini_section, 'candidate_genes', required=False) else "",
            output=output
        )
    )

def html(
    input_tsv,
    output_dir,
    output_prefix,
    ini_section='annotSV'
    ):
    """
    Create html output from AnnotSV tsv.

    :return: a job for AnnotSV html creation
    """
    output = [
        os.path.join(output_dir, output_prefix + ".html")
    ]

    return Job(
        [input_tsv],
        output,
        [
            [ini_section, "module_perl"],
            [ini_section, "module_knotAnnotSV"]
        ],
        command="""\
perl $KNOTANNOTSV_HOME/knotAnnotSV.pl {html_options} \\
  --configFile {config} \\
  --datatableDir $KNOTANNOTSV_HOME/DataTables \\
  --annotSVfile {input_tsv} \\
  -genomeBuild {genome_build} \\
  --outDir {output_dir} \\
  --outPrefix {output_prefix}""".format(
            html_options=global_conf.global_get(ini_section, 'html_options', required=False),
            config=global_conf.global_get(ini_section, 'knotAnnotSV_config', required=False) if global_conf.global_get(ini_section, 'knotAnnotSV_config', required=False) else "$KNOTANNOTSV_CONFIG",
            input_tsv=input_tsv,
            genome_build=global_conf.global_get(ini_section, 'assembly'),
            output_dir=output_dir,
            output_prefix=output_prefix
        )
    )

def excel(
    input_tsv,
    output_dir,
    output_prefix,
    ini_section='annotSV'
    ):
    """
    Create excel output from AnnotSV tsv.

    :return: a job for AnnotSV html creation
    """
    output = [
        os.path.join(output_dir, output_prefix + ".xlsm")
    ]

    return Job(
        [input_tsv],
        output,
        [
            [ini_section, "module_perl"],
            [ini_section, "module_knotAnnotSV"]
        ],
        command="""\
perl $KNOTANNOTSV_HOME/knotAnnotSV2XL.pl {excel_options} \\
  --configFile {config} \\
  --datatableDir $KNOTANNOTSV_HOME/DataTables \\
  --annotSVfile {input_tsv} \\
  -genomeBuild {genome_build} \\
  --geneCountThreshold {genecount_threshold} \\
  --outDir {output_dir} \\
  --outPrefix {output_prefix}""".format(
            excel_options=global_conf.global_get(ini_section, 'excel_options', required=False),
            config=global_conf.global_get(ini_section, 'knotAnnotSV_config', required=False) if global_conf.global_get(ini_section, 'knotAnnotSV_config', required=False) else "$KNOTANNOTSV_CONFIG",
            input_tsv=input_tsv,
            genome_build=global_conf.global_get(ini_section, 'assembly'),
            genecount_threshold=global_conf.global_get(ini_section, 'genecount_threshold'),
            output_dir=output_dir,
            output_prefix=output_prefix
        )
    )