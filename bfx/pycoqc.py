#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *


def pycoqc(readset_name,
           input_summary,
           output_directory,
           input_bam
           ):
    """
    Create a pycoQC job for nanopore reads QC.

    :return: a job for nanopore QC using pycoQC
    """

    min_qual = config.param('pycoqc', 'min_pass_qual')

    in_bam = ' --bam_file ' + input_bam

    out_html = readset_name + ".html"
    out_json = readset_name + ".json"

    return Job(
        [input_summary],
        [os.path.join(output_directory, out_html),
         os.path.join(output_directory, out_json)],
        [["pycoqc", "module_python3"]],
        command="""\
mkdir -p {output_directory} && \\
cd {output_directory} && \\
pycoQC --verbose {other_options} --summary_file {input_summary}{in_bam} \\
    --report_title {readset_name} \\
    --min_pass_qual {min_qual} \\
    --html_outfile {out_html} \\
    --json_outfile {out_json}
        """.format(
            output_directory=output_directory,
            other_options=config.param('pycoqc', 'other_options', required=False),
            input_summary=input_summary,
            in_bam=in_bam,
            readset_name=readset_name,
            min_qual=min_qual,
            out_html=out_html,
            out_json=out_json
        )
    )
