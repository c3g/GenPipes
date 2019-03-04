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
import logging
import os

# MUGQIC Modules
from core.config import *
from core.job import * 

def bamqc(input_bam, output_directory, output, options):

    inputs = [input_bam]
    outputs = [output]

    (input_basename, file_format) = os.path.splitext(input_bam)

    return Job(
        inputs,
        outputs,
        [
            ['qualimap_bamqc', 'module_java'],
            ['qualimap_bamqc', 'module_qualimap'],
        ],
        command="""\
qualimap bamqc {other_options} \\
  -bam {input_bam} -outdir {output_directory} \\
  --java-mem-size={ram}""".format(
            input_bam=input_bam,
            output_directory=output_directory,
            other_options=options,
            ram=config.param('qualimap_bamqc', 'ram'),
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
            ['qualimap_rnaseq', 'module_java'],
            ['qualimap_rnaseq', 'module_qualimap'],
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
            gtf=config.param('qualimap', 'gtf', type='filepath'),
            output_directory=output_directory,
            output=output,
            ram=config.param('qualimap', 'ram'),
            other_options=config.param('qualimap_rnaseq', 'other_options')
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

#Job(command="""\
#unset DISPLAY \\
#R --no-save --no-restore<<-'EOF'
  #fns=list.dirs("{input_directory}",full.names=T,recursive=F)
  #write.table(data.frame(basename(fns), file.path(fns,"bamqc")), file="{output_directory}/multi-bamqc_list.txt", sep='\t', quote=F, col.names=F, row.names=F)
#EOF""".format(
                #input_directory=input_directory,
                #output_directory=output_directory
            #)),

    return job
