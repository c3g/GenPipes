#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def exploratory_rnaseq(input_file, genes_file, output_dir):

    return Job(
        [genes_file],
        [],
        [['gq_seq_utils_exploratory_rnaseq', 'module_R']],
        command = """\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(gqSeqUtils))

exploratoryRNAseq(input.path="{input_file}", genes.path="{genes_file}", output.dir="{output_dir}")
print("done.")

EOF""".format(
        input_file=input_file,
        genes_file=genes_file,
        output_dir=output_dir
    ))

def report(ini_filepaths, project_path, pipeline_type, output_directory):

    title = config.param('gq_seq_utils_report', 'report_title', required=False)
    path = os.path.join(output_directory, config.param('gq_seq_utils_report', 'report_dir'))
    author = config.param('gq_seq_utils_report', 'report_author', required=False)
    contact = config.param('gq_seq_utils_report', 'report_contact', required=False)

    # Job input files must be set in pipeline class since they are different for each pipeline
    job = Job([], [os.path.join(path, "index.html")], [['gq_seq_utils_report', 'module_R']])

    job.command = """\
R --no-save -e 'library(gqSeqUtils); mugqicPipelineReport(pipeline="{pipeline}"{title}{path}{author}{contact}, ini.file.path=c({ini_filepaths}), project.path="{project_path}")'""".format(
        pipeline=pipeline_type,
        title=", report.title=\"" + title + "\"" if title else "",
        path=", report.path=\"" + path + "\"" if path else "",
        author=", report.author=\"" + author + "\"" if author else "",
        contact=", report.contact=\"" + contact + "\"" if contact else "",
        ini_filepaths=",".join(['"' + ini_filepath + '"' for ini_filepath in ini_filepaths]),
        project_path=project_path
    )

    return job
