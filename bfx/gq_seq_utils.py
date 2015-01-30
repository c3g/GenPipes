#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def exploratory_analysis_rnaseq(htseq_count_file, cuffnorm_dir, genes_file, output_dir):

    return Job(
        [htseq_count_file, os.path.join(cuffnorm_dir, "isoforms.fpkm_table"), os.path.join(cuffnorm_dir, "isoforms.attr_table")],
        [os.path.join(output_dir, "top_sd_heatmap_cufflinks_logFPKMs.pdf")],
        [
            ['gq_seq_utils_exploratory_analysis_rnaseq', 'module_R'],
            ['gq_seq_utils_exploratory_analysis_rnaseq', 'module_mugqic_R_packages']
        ],
        command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(gqSeqUtils))

exploratoryAnalysisRNAseq(htseq.counts.path="{htseq_count_file}", cuffnorm.fpkms.dir="{cuffnorm_dir}", genes.path="{genes_file}", output.dir="{output_dir}")
print("done.")

EOF""".format(
        htseq_count_file=htseq_count_file,
        cuffnorm_dir=cuffnorm_dir,
        genes_file=genes_file,
        output_dir=output_dir
    ))

def report(ini_filepaths, project_path, pipeline_type, output_directory):

    title = config.param('gq_seq_utils_report', 'report_title', required=False)
    path = os.path.join(output_directory, config.param('gq_seq_utils_report', 'report_dir'))
    author = config.param('gq_seq_utils_report', 'report_author', required=False)
    contact = config.param('gq_seq_utils_report', 'report_contact', required=False)

    return Job(
        # Job input files must be set in pipeline class since they are different for each pipeline
        [],
        [os.path.join(path, "index.html")],
        [
            ['gq_seq_utils_report', 'module_R'],
            ['gq_seq_utils_report', 'module_mugqic_R_packages']
        ],
        command="""\
R --no-save -e 'library(gqSeqUtils); mugqicPipelineReport(pipeline="{pipeline}"{title}{path}{author}{contact}, ini.file.path=c({ini_filepaths}), project.path="{project_path}")'""".format(
        pipeline=pipeline_type,
        title=", report.title=\"" + title + "\"" if title else "",
        path=", report.path=\"" + path + "\"" if path else "",
        author=", report.author=\"" + author + "\"" if author else "",
        contact=", report.contact=\"" + contact + "\"" if contact else "",
        ini_filepaths=",".join(['"' + os.path.abspath(ini_filepath) + '"' for ini_filepath in ini_filepaths]),
        project_path=project_path
    ))
