#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def blastn_on_db(db, query, output, other_options=""):
    job = Job([query], [output], [['DEFAULT', 'module_blast']])

    job.command = """\
blastn {other_options} \\
  -db {db} \\
  -query {query} \\
  -out {output}""".format(
        other_options=other_options,
        db=db,
        query=query,
        output=output
    )

    return job
