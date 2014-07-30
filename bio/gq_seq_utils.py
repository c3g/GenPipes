#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def report(ini_filepath, project_path, pipeline_type, output_directory):

    title = config.param('report', 'project_name', required=False)
    path = os.path.join(output_directory, config.param('report', 'report_name'))
    author = config.param('report', 'report_author', required=False)
    contact = config.param('report', 'report_contact', required=False)

    # Job input files must be set in pipeline since they are different in each pipeline
    job = Job([], [os.path.join(path, "index.html")], [['report', 'module_R']])

    job.command = \
"""R --no-save -e 'library(gqSeqUtils); mugqicPipelineReport(pipeline=\\"{pipeline}\\"{title}{path}{author}{contact}, ini.file.path=\\"{ini_filepath}\\", project.path=\\"{project_path}\\")'""".format(
        pipeline=pipeline_type,
        title=", report.title=\\\"" + title + "\\\"" if title else "",
        path=", report.path=\\\"" + path + "\\\"" if path else "",
        author=", report.author=\\\"" + author + "\\\"" if author else "",
        contact=", report.contact=\\\"" + contact + "\\\"" if contact else "",
        ini_filepath=ini_filepath,
        project_path=project_path
    )

    return job
