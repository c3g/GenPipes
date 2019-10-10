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
import re

# MUGQIC Modules
from core.config import config
from core.job import Job

def index(
    input,
    output,
    basecalls_dir,
    mismatches,
    lane,
    mask
    ):

    return Job(
        [input],
        [output],
        [
            ["index", "module_java"],
            ["index", "module_mugqic_tools"]
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} \\
  {java_other_options} \\
  -Xmx{ram} \\
  -jar $JAVA_TOOLS/{jar} \\
  MAX_MISMATCHES={mismatches} \\
  NUM_PROCESSORS={threads} \\
  BARCODE_FILE={barcode_file} \\
  BASECALLS_DIR={basecalls_dir} \\
  LANE={lane_number} \\
  READ_STRUCTURE={read_structure} \\
  METRICS_FILE={output} \\
  TMP_DIR={tmp_dir}""".format(
            tmp_dir=config.param('index', 'tmp_dir'),
            java_other_options=config.param('index', 'java_other_options'),
            ram=config.param('index', 'ram'),
            jar=config.param('index', 'jar'),
            mismatches=mismatches,
            threads=config.param('index', 'threads'),
            barcode_file=config.param('index', 'barcode_file'),
            basecalls_dir=basecalls_dir,
            lane_number=lane,
            read_structure=mask,
            output=output
        )
    )

def bcl2fastq(
    input,
    fastq_outputs,
    output_dir,
    sample_sheet,
    demultiplexing,
    run,
    lane,
    mismatches,
    mask
    ):

    command_suffix = ""
    if demultiplexing:
        command_suffix = """\
  --barcode-mismatches {number_of_mismatches} \\
  --use-bases-mask {mask}""".format(
            number_of_mismatches=mismatches,
            mask=mask
        )

    return Job(
        [input],
        fastq_outputs,
        [
            ['fastq', 'module_bcl_to_fastq']
        ],
        command="""\
bcl2fastq \\
  --runfolder-dir {run_dir} \\
  --output-dir {output_dir} \\
  --tiles {tiles} \\
  --sample-sheet {sample_sheet} \\
  --create-fastq-for-index-reads \\
  {additional_parameters} \\
  {other_options}""".format(
            run_dir=run,
            output_dir=output_dir,
            tiles="s_" + str(lane),
            sample_sheet=sample_sheet,
            additional_parameters=command_suffix,
            other_options=config.param('fastq', 'other_options')
        ),
        name="fastq." + self.run_id + "." + str(self.lane_number),
        samples=self.samples
    )
