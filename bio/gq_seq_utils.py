#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def client_report(ini_filepath, project_path, pipeline_type):
    # Job input and output files must be set in pipeline since they are different in each pipeline
    job = Job([], [], [['report', 'moduleVersion.R']])

    title = config.param('report', 'projectName', required=False)
    path = config.param('report', 'report.path', required=False)
    author = config.param('report', 'report.author', required=False)
    contact = config.param('report', 'report.contact', required=False)

    job.command = \
"""R --no-save -e 'library(gqSeqUtils); mugqicPipelineReport( \\
  pipeline=\\"{pipeline}\\"{title}{path}{author}{contact}, \\
  ini.file.path=\\"{ini_filepath}\\", \\
  project.path=\\"{project_path}\\" \\
)'""".format(
        pipeline=pipeline_type,
        title=", \\\n  report.title=\\\"" + title + "\\\"" if title else "",
        path=", \\\n  report.path=\\\"" + path + "\\\"" if path else "",
        author=", \\\n  report.author=\\\"" + author + "\\\"" if author else "",
        contact=", \\\n  report.contact=\\\"" + contact + "\\\"" if contact else "",
        ini_filepath=ini_filepath,
        project_path=project_path
    )

    return job
