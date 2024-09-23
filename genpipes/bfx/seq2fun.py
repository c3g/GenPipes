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


def processing(input_files, output_file, sample_file, profiling, ini_section='seq2fun'):
    genemap = global_conf.global_get(ini_section, 'genemap')
    tfmi = global_conf.global_get(ini_section, 'tfmi')
    other_options = global_conf.global_get(ini_section, 'other_options')
    return Job(
        input_files,
        output_file,
        [
            [ini_section, 'module_seq2fun']
        ],
        command="""seq2fun --sampletable {sample_file} --tfmi {tfmi} --genemap {genemap} {profiling} {other_options}""".format(
            sample_file=sample_file,
            tfmi = tfmi,
            genemap = genemap,
            other_options = other_options,
            profiling = profiling
                   )

    )

def deseq2(
    design_file,
    count_matrix,
    output_dir, 
    ini_section='seq2fun'
    ):

    localfit = "-l" if global_conf.global_get('differential_expression_deseq', 'localfit') else ""

    return  Job(
        [count_matrix],
        [os.path.join(output_dir, "deseq_results.csv"), os.path.join(output_dir, "dge_results.csv")],
        [
            [ini_section, 'module_R'],
            [ini_section, 'module_mugqic_tools'],
        ],
        command="""\
Rscript $R_TOOLS/deseq2.R \\
  -d {design_file} \\
  -c {count_matrix} \\
  -o {output_dir} \\
  {localfit}""".format(
        design_file=design_file,
        count_matrix=count_matrix,
        output_dir=output_dir,
        localfit=localfit
    ))


def ko_pathway_analysis(diff_report, output_prefix, output_dir, ini_section='seq2fun_pathway'):
    fdr = global_conf.global_get(ini_section, 'fdr')
    rds_file = global_conf.global_get(ini_section, 'rds')
    map_list = global_conf.global_get(ini_section, 'user_pathway_list')
    kegg_all = global_conf.global_get(ini_section, 'kegg_all')
    return Job(
        [diff_report],
        [os.path.join(output_dir, output_prefix + ".txt")],
        [
            ['seq2fun', 'module_R'],
            ['seq2fun', 'module_mugqic_tools']
        ],
        command="""\
Rscript $R_TOOLS/KOPathawayAnalysis.R \\
    -i {diff_report} \\
    -map {map_list} \\
    -o {output_dir} \\
    -p {output_prefix} \\
    -rds {rds_file} \\
    -kegg {kegg_all} \\
    -fdr {fdr}""".format(
            diff_report=diff_report,
            map_list=map_list,
            output_dir=output_dir,
            fdr=fdr,
            output_prefix=output_prefix,
            rds_file=rds_file,
            kegg_all=kegg_all
        ))