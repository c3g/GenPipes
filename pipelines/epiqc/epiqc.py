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
from core.job import *
from core.config import *

def modify_inputinfofile(input_file, sample, histone,  output_inputinfo_file):
    return Job(
        [input_file, output_inputinfo_file],
        [output_inputinfo_file],
        [],
        command="""\
        echo -e "{sample}\\t{histone}\\t{sample}.bedgraph.gz" >> {inputinfofile}
        """.format(
            input_files=input_file,
            sample=sample,
            histone=histone,
            inputinfofile=output_inputinfo_file

        )
    )


def convert(input_dir, output_dir, output_files, inputinfofile, histone_mark, sample):
    # input = input_files.extend(inputinfofile)
    return Job(
        [inputinfofile],
        output_files,
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name="chromimpute_convert." + sample + "." + histone_mark,
        command="""\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
  Convert {chrom} \\
  -m {histone_mark} \\
  -l {convertsample} \\
  -r {resolution} \\
  {path_to_dataset} \\
  {inputinfofile} \\
  {chrom_sizes} \\
  {output_dir}""".format(
      java_other_options=config.param('DEFAULT', 'java_other_options'),
      ram=config.param('chromimpute', 'ram'),
      chrom="-c " + config.param('chromimpute', 'chrom') if config.param('chromimpute', 'chrom') else "",
      histone_mark=histone_mark,
      convertsample=sample,
      resolution=config.param('chromimpute', 'resolution'),
      path_to_dataset=input_dir,
      inputinfofile=inputinfofile,
      chrom_sizes="<(awk '{print $1\"\\t\"$2}' %s)" % (config.param('chromimpute', 'chromosome_size')),# config.param('chromimpute', 'chrominfofile')
      output_dir=output_dir
      )
    )

def compute_global_dist(input_files, output_dir, output_files, converteddir, inputinfofile, histone_mark):
    return Job(
        input_files,
        output_files,
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name="chromimpute_compute_global_dist." + histone_mark,
        command="""\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
  ComputeGlobalDist \\
  -m {histone_mark} \\
  -r {resolution} \\
  {converteddir} \\
  {inputinfofile} \\
  {chrom_sizes} \\
  {output_dir}""".format(
      java_other_options=config.param('DEFAULT', 'java_other_options'),
      ram=config.param('chromimpute', 'ram'),
      histone_mark=histone_mark,
      resolution=config.param('chromimpute', 'resolution'),
      converteddir=converteddir,
      inputinfofile=inputinfofile,
      chrom_sizes="<(awk '{print $1\"\\t\"$2}' %s)" % (config.param('chromimpute', 'chromosome_size')),# config.param('chromimpute', 'chrominfofile')
      output_dir=output_dir
      )
    )

def generate_train_data(input_files, output_dir, output_files, converteddir, distancedir, inputinfofile, histone_mark):
    return Job(
        input_files,
        output_files,
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name="chromimpute_generate_train_data." + histone_mark,
        command="""\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
  GenerateTrainData \\
  -r {resolution} \\
  {converteddir} \\
  {distancedir} \\
  {inputinfofile} \\
  {chrom_sizes} \\
  {output_dir} \\
  {histone_mark}""".format(
      java_other_options=config.param('DEFAULT', 'java_other_options'),
      ram=config.param('chromimpute', 'ram'),
      # chrom=chrom,#ddddd.param('chromimpute', 'chrom'),
      resolution=config.param('chromimpute', 'resolution'),
      converteddir=converteddir,
      distancedir=distancedir,
      inputinfofile=inputinfofile,
      chrom_sizes="<(awk '{print $1\"\\t\"$2}' %s)" % (config.param('chromimpute', 'chromosome_size')),# config.param('chromimpute', 'chrominfofile')
      output_dir=output_dir,
      histone_mark=histone_mark
      )
    )


def train(input_files, output_dir, output_files, traindatadir, inputinfofile, sample, histone_mark):
    return Job(
        input_files,
        output_files,
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name="chromimpute_train." + sample + "_" + histone_mark,
        command="""\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
  Train \\
  {traindatadir} \\
  {inputinfofile} \\
  {predictordir} \\
  {sample} \\
  {histone_mark}""".format(
      java_other_options=config.param('DEFAULT', 'java_other_options'),
      ram=config.param('chromimpute', 'ram'),
      traindatadir=traindatadir,
      inputinfofile=inputinfofile,
      predictordir=output_dir,
      sample=sample,
      histone_mark=histone_mark
      )
    )

def apply(input_files, output_dir, converteddir, distancedir, predictordir, inputinfofile, sample, mark):
    return Job(
        input_files,
        [output_dir],
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name="chromimpute_apply."+sample+"_"+mark,
        command="""\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
    Apply \\
    -r {resolution} \\
    {converteddir} \\
    {distancedir} \\
    {predictordir} \\
    {inputinfofile} \\
    {chrom_sizes} \\
    {output_dir} \\
    {sample} \\
    {mark}""".format(
        java_other_options=config.param('DEFAULT', 'java_other_options'),
        ram=config.param('chromimpute', 'ram'),
        # chrom=config.param('chromimpute', 'chrom'),
        resolution=config.param('chromimpute', 'resolution'),
        converteddir=converteddir,
        distancedir=distancedir,
        predictordir=predictordir,
        inputinfofile=inputinfofile,
        chrom_sizes="<(awk '{print $1\"\\t\"$2}' %s)" % (config.param('chromimpute', 'chromosome_size')),# config.param('chromimpute', 'chrominfofile')
        output_dir=output_dir,
        sample=sample,
        mark=mark
        )
    )

def eval(input_dir, percent1, percent2, converteddir, converted_file, output_dir, imputed_file, output_path):
    return Job(
        [input_dir],
        [output_path],
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name="chromimpute_eval."+converted_file+"."+imputed_file,
        command="""\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
    Eval \\
    -p {percent1} {percent2} \\
    {converteddir} \\
    {converted_file} \\
    {output_dir} \\
    {imputed_file} \\
    {chrom_sizes} > {output_path}""".format(
        java_other_options=config.param('DEFAULT', 'java_other_options'),
        ram=config.param('chromimpute', 'ram'),
        percent1=percent1,
        percent2=percent2,
        converteddir=converteddir,
        converted_file=converted_file,
        output_dir=output_dir,
        imputed_file=imputed_file,
        chrom_sizes=config.param('chromimpute', 'chrominfofile'),
        output_path=output_path
        )
    )
