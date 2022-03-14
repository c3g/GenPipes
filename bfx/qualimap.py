################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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
import logging
import os

# MUGQIC Modules
from core.config import *
from core.job import * 

def bamqc(input_bam, output_directory, output, options, ini_section='qualimap'):

    inputs = [input_bam]
    outputs = [output]

    (input_basename, file_format) = os.path.splitext(input_bam)

    return Job(
        inputs,
        outputs,
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_qualimap'],
        ],
        command="""\
qualimap bamqc {other_options} \\
  -bam {input_bam} -outdir {output_directory} \\
  --java-mem-size={ram}""".format(
            input_bam=input_bam,
            output_directory=output_directory,
            other_options=options,
#            bed="\\\n  --feature-file " + bed if bed else "",
            ram=global_config_parser.param(ini_section, 'ram'),
        ),
        removable_files=[]
    )

def rnaseq(input_bam, output_directory, output):

    inputs = [input_bam]

    outputs = [output]

    (input_basename, file_format) = os.path.splitext(input_bam)

    return Job(
        inputs,
        outputs,
        [
            ['DEFAULT', 'module_java'],
            ['qualimap', 'module_qualimap'],
        ],
        command="""\
qualimap rnaseq \\
  -bam {input_bam} \\
  -gtf {gtf} \\
  -outdir {output_directory} \\
  -oc {output} \\
  --java-mem-size={ram} \\
  {other_options}""".format(
            input_bam=input_bam,
            gtf=global_config_parser.param('qualimap', 'gtf', param_type='filepath'),
            output_directory=output_directory,
            output=output,
            ram=global_config_parser.param('qualimap', 'ram'),
            other_options=global_config_parser.param('qualimap', 'other_options')
        ),
        removable_files=[]
    )

def multibamqc(inputs, output_directory):

    outputs = [os.path.join(output_directory, "report.html")]

    job=Job(
        inputs,
        outputs,
        [
            ['qualimap_multibamqc', 'module_qualimap'],
            ['qualimap_multibamqc', 'module_R']
        ],
        command="""\
qualimap multi-bamqc \\
  -d {output_directory}/multi-bamqc_list.txt \\
  -outdir {output_directory} \\
  -outfile {outfile}""".format(
            output_directory=output_directory,
            outfile=os.path.join(output_directory, "report.html")
        ),
        removable_files=[]
    )

    job=concat_jobs([
        Job(command="""\
for i in {input_files}; do \\
  path1=$(dirname $i); \\
  path2=$(dirname $path1); \\
  echo -e \"$(basename $path2)\t$path1\"; \\
done > {output_directory}/multi-bamqc_list.txt""".format(
                input_files=" ".join(inputs),
                output_directory=output_directory
            )),
        job
    ])

    return job
