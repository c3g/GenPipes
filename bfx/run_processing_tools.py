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
import sys
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

    barcode_file = config.param('index', 'barcode_file', type='filepath', required='false')
    if not (barcode_file and os.path.isfile(barcode_file)):
        barcode_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'resources', 'barcode.mergedup.txt')

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
            barcode_file=barcode_file,
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
    run,
    lane,
    extra_option,
    demultiplex=False,
    mismatches=None,
    mask=None,
    ):

    if demultiplex:
        demultiplex_parameters = """\
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
  {demultiplex_parameters} \\
  {other_options} \\
  {extra_option}""".format(
            run_dir=run,
            output_dir=output_dir,
            tiles="s_" + str(lane),
            sample_sheet=sample_sheet,
            demultiplex_parameters=demultiplex_parameters if demultiplex_parameters else "",
            other_options=config.param('fastq', 'other_options'),
            extra_option=extra_option
        )
    )

def fastqmultx(
    barcode_file,
    mismatches,
    outputs,
    R1_fastq,
    I1_fastq,
    R2_fastq,
    I2_fastq=None
    ):

    metrics_file = outputs[0]
    output_dir = os.path.dirname(metrics_file)
    return Job(
        [
            R1_fastq,
            I1_fastq,
            R2_fastq,
            I2_fastq if I2_fastq else None
        ],
        outputs,
        [
            ['fastq', 'module_fastq_multx']
        ],
        command="""\
fastq-multx \\
  -B {barcode} \\
  -m {mismatches} \\
  {input_I1} \\
  {input_I2} \\
  {input_R1} \\
  {input_R2} \\
  {output_I1_fake} \\
  {output_I2_fake} \\
  {output_R1} \\
  {output_R2} \\
  > {stdout_metrics}""".format(
            barcode=barcode_file,
            mismatches=mismatches,
            input_I1=I1_fastq,
            input_I2=I2_fastq if I2_fastq else "",
            input_R1=R1_fastq,
            input_R2=R2_fastq,
            output_I1_fake="-o " + os.path.join(output_dir, "%_I1.fastq"),
            output_I2_fake="-o " + os.path.join(output_dir, "%_I2.fastq") if I2_fastq else "",
            output_R1="-o " + os.path.join(output_dir, "%_R1.fastq"),
            output_R2="-o " + os.path.join(output_dir, "%_R2.fastq"),
            stdout_metrics=metrics_file
        ),
        removable_files=outputs+[output_dir]
    )

def demux_fastqs(
    sample_sheet,
    mismatches,
    mask,
    outputs,
    R1_fastq,
    R2_fastq
    ):

    output_dir = os.path.dirname(outputs[0])
    metrics_file = os.path.join(output_dir, "DemuxFastqs.metrics.txt")
    return Job(
        [
            R1_fastq,
            R2_fastq
        ],
        outputs + [ metrics_file],
        [
            ['fastq', 'module_java'],
            ['fastq', 'module_fgbio']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} \\
  {java_other_options} \\
  -Xmx{ram} \\
  -jar $FGBIO_JAR DemuxFastqs \\
  --threads {threads} \\
  --max-mismatches {mismatches} \\
  --metrics {metrics_file} \\
  --inputs {inputs} \\
  --read-structures {read_structure} \\
  --metadata {sample_sheet} \\
  --output {output_dir} \\
  --output-type fastq \\
  --include-all-bases-in-fastqs true""".format(
            tmp_dir=config.param('fastq', 'tmp_dir'),
            java_other_options=config.param('fastq', 'java_other_options'),
            ram=config.param('fastq', 'ram'),
            threads=config.param('fastq', 'threads'),
            mismatches=mismatches,
            metrics_file=metrics_file,
            inputs=" ".join([R1_fastq, R2_fastq]),
            read_structure=mask,
            sample_sheet=sample_sheet,
            output_dir=output_dir
        )
    )

