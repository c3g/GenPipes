#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def deseq(
    design_file,
    count_matrix,
    output_dir
    ):

    return  Job(
        [count_matrix],
        [output_dir],
        [
            ['differential_expression_deseq', 'module_mugqic_tools'],
            ['differential_expression_deseq', 'module_R']
        ],
        command = """\
Rscript $R_TOOLS/deseq.R \\
  -d {design_file} \\
  -c {count_matrix} \\
  -o {output_dir}""".format(
        design_file=design_file,
        count_matrix=count_matrix,
        output_dir=output_dir
    ))

def edger(
    design_file,
    count_matrix,
    output_dir
    ):

    return  Job(
        [count_matrix],
        [output_dir],
        [
            ['differential_expression_edger', 'module_mugqic_tools'],
            ['differential_expression_edger', 'module_R']
        ],
        command = """\
Rscript $R_TOOLS/edger.R \\
  -d {design_file} \\
  -c {count_matrix} \\
  -o {output_dir}""".format(
        design_file=design_file,
        count_matrix=count_matrix,
        output_dir=output_dir
    ))

def goseq(
    result_file,
    output_file,
    columns,
    method
    ):

    max_go_results = config.param('differential_expression_goseq', 'max_go_results', required=False, type='posint')
    gene_size_file = config.param('differential_expression_goseq', 'gene_size', required=False, type='filepath')
    go_link_file = config.param('differential_expression_goseq', 'go_link_file', required=False, type='filepath')
    gene_id_type = config.param('differential_expression_goseq', 'gene_id_type', required=False)

    return  Job(
        [result_file],
        [output_file],
        [
            ['differential_expression_goseq', 'module_mugqic_tools'],
            ['differential_expression_goseq', 'module_R']
        ],
        command = """\
Rscript $R_TOOLS/goseq.R {max_go_results}{gene_size_file}{go_link_file}{gene_id_type} \\
  -d {result_file} \\
  -c {columns} \\
  -t {go_annotation} \\
  -k {geneid2symbol} \\
  -s {goseq_genome} \\
  -o {output_file}""".format(
        max_go_results="\\\n  -m " + max_go_results if max_go_results else "",
        gene_size_file="\\\n  -a " + gene_size_file if gene_size_file else "",
        go_link_file="\\\n  -G " + go_link_file if go_link_file else "",
        gene_id_type="\\\n  -i " + gene_id_type if gene_id_type else "",
        result_file=result_file,
        columns=columns,
        output_file=output_file
    ))
