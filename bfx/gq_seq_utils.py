#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def exploratory_rnaseq(readset_filepath, project_path, ini_filepaths):

    job = Job([], [], [['gq_seq_utils_exploratory_rnaseq', 'module_R']])

    job.command = """\
R --no-save <<EOF
suppressPackageStartupMessages(library(gqSeqUtils))

initIllmSeqProject(nanuq.file="{readset_filepath}", overwrite.sheets=TRUE, project.path="{project_path}")
exploratoryRNAseq(project.path="{project_path}", ini.file.path=c({ini_filepaths}))
print("done.")

EOF""".format(
        readset_filepath=readset_filepath,
        project_path=project_path,
        ini_filepaths=",".join(['"' + ini_filepath + '"' for ini_filepath in ini_filepaths])
    )

    return job

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
