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

import os

from core.job import *
from core.config import *

def convert(input_dir, output_dir, inputinfofile, mark):
    return Job(
        [input_dir],
        [output_dir],
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name = "chromimpute_convert_"+mark,
        command = """\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
  Convert \\
  {chrom} \\
  -m {mark} \\
  {resolution} \\
  {path_to_dataset} \\
  {inputinfofile} \\
  {chrom_sizes} \\
  {output_dir}""".format(
        java_other_options = config.param('DEFAULT', 'java_other_options'),
        ram = config.param('chromimpute', 'ram'),
        chrom = config.param('chromimpute','chrom'),
        mark = mark,
        resolution = config.param('chromimpute', 'resolution'),
        path_to_dataset = config.param('chromimpute', 'dataset'),
        inputinfofile = inputinfofile,
        chrom_sizes = config.param('chromimpute', 'chromsizes'),
        output_dir = output_dir
        )
    )

def compute_global_dist(input_dir, output_dir, inputinfofile, mark):
    return Job(
        [input_dir],
        [output_dir],
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name = "chromimpute_compute_global_dist_"+mark,
        command = """\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
  ComputeGlobalDist \\
  -m {mark} \\
  {resolution} \\
  {converteddir} \\
  {inputinfofile} \\
  {chrom_sizes} \\
  {output_dir}""".format(
        java_other_options = config.param('DEFAULT', 'java_other_options'),
        ram = config.param('chromimpute', 'ram'),
        mark = mark,
        resolution = config.param('chromimpute', 'resolution'),
        converteddir = input_dir,
        inputinfofile = inputinfofile,
        chrom_sizes = config.param('chromimpute', 'chromsizes'),
        output_dir = output_dir
        )
    )

def generate_train_data(input_dir, output_dir, converteddir, inputinfofile, mark):
    return Job(
        [input_dir],
        [output_dir],
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name = "chromimpute_generate_train_data_"+mark,
        command = """\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
  GenerateTrainData \\
  {chrom} \\
  {resolution} \\
  {converteddir} \\
  {distancedir} \\
  {inputinfofile} \\
  {chrom_sizes} \\
  {output_dir} \\
  {mark}""".format(
        java_other_options = config.param('DEFAULT', 'java_other_options'),
        ram = config.param('chromimpute', 'ram'),
        chrom = config.param('chromimpute','chrom'),
        resolution = config.param('chromimpute', 'resolution'),
        converteddir = converteddir,
        distancedir = input_dir,
        inputinfofile = inputinfofile,
        chrom_sizes = config.param('chromimpute', 'chromsizes'),
        output_dir = output_dir,
        mark = mark
        )
    )


def train(input_dir, output_dir, inputinfofile, sample, mark):
    return Job(
        [input_dir],
        [output_dir],
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name = "chromimpute_train_"+sample+"_"+mark,
        command = """\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
  Train \\
  {traindatadir} \\
  {inputinfofile} \\
  {predictordir} \\
  {sample} \\
  {mark}""".format(
        java_other_options = config.param('DEFAULT', 'java_other_options'),
        ram = config.param('chromimpute', 'ram'),
        traindatadir = input_dir,
        inputinfofile = inputinfofile,
        predictordir = output_dir,
        sample = sample,
        mark = mark
        )
    )

def apply(input_dir, output_dir, converteddir, distancedir, predictordir, inputinfofile, sample, mark):
    return Job(
        [input_dir],
        [output_dir],
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name = "chromimpute_apply_"+sample+"_"+mark,
        command = """\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
    Apply \\
    {chrom} \\
    {resolution} \\
    {converteddir} \\
    {distancedir} \\
    {predictordir} \\
    {inputinfofile} \\
    {chrom_sizes} \\
    {output_dir} \\
    {sample} \\
    {mark}""".format(
        java_other_options = config.param('DEFAULT', 'java_other_options'),
        ram = config.param('chromimpute', 'ram'),
        chrom = config.param('chromimpute','chrom'),
        resolution = config.param('chromimpute', 'resolution'),
        converteddir = converteddir,
        distancedir = distancedir,
        predictordir = predictordir,
        inputinfofile = inputinfofile,
        chrom_sizes = config.param('chromimpute', 'chromsizes'),
        output_dir = output_dir,
        sample = sample,
        mark = mark
        )
    )

def eval(input_dir, percent1, percent2, converteddir, converted_file, imputed_file, output_path):
    return Job(
        [input_dir],
        [output_path],
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name = "chromimpute_eval_"+converted_file+"_vs_"+imputed_file,
        command = """\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
    Eval \\
    -p {percent1} {percent2} \\
    {converteddir} \\
    {converted_file} \\
    {input_dir} \\
    {imputed_file} \\
    {chrom_sizes} > {output_path}""".format(
        java_other_options = config.param('DEFAULT', 'java_other_options'),
        ram = config.param('chromimpute', 'ram'),
        percent1 = percent1,
        percent2 = percent2,
        converteddir = converteddir,
        converted_file = converted_file,
        input_dir = input_dir,
        imputed_file = imputed_file,
        chrom_sizes = config.param('chromimpute', 'chromsizes'),
        output_path = output_path
        )
    )



