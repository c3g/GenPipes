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

def convert(input_dir, output_dir):
    return Job(
        [input_dir],
        [output_dir],
        [['java', 'module_java']],
        name = "chromimpute_convert",
        command = """\
java {java_options} -jar /lb/project/mugqic/projects/rami_test/Tools/chromimpute/ChromImpute.jar \\
  Convert \\
  -c {chrom} \\
  -r {resolution} \\
  {path_to_dataset} \\
  {inputinfofile} \\
  {chrom_sizes} \\
  {output_dir}""".format(
        java_options = config.param('DEFAULT', 'java_options'),
        chrom = config.param('chromimpute','chrom'),
        resolution = config.param('chromimpute', 'resolution'),
        path_to_dataset = config.param('chromimpute', 'dataset'),
        inputinfofile = config.param('chromimpute', 'inputinfofile'),
        chrom_sizes = config.param('chromimpute', 'chromsizes'),
        output_dir = output_dir
        )
    )

def compute_global_dist(input_dir, output_dir):
    return Job(
        [input_dir],
        [output_dir],
        [['java', 'module_java']],
        name = "chromimpute_compute_global_dist",
        command = """\
java {java_options} -jar /lb/project/mugqic/projects/rami_test/Tools/chromimpute/ChromImpute.jar \\
  ComputeGlobalDist \\
  -r {resolution} \\
  {converteddir} \\
  {inputinfofile} \\
  {chrom_sizes} \\
  {output_dir}""".format(
        java_options = config.param('DEFAULT', 'java_options'),
        resolution = config.param('chromimpute', 'resolution'),
        converteddir = input_dir,
        inputinfofile = config.param('chromimpute', 'inputinfofile'),
        chrom_sizes = config.param('chromimpute', 'chromsizes'),
        output_dir = output_dir
        )
    )


def generate_train_data(input_dir, converteddir, output_dir, mark):
    return Job(
        [input_dir],
        [output_dir],
        [['java', 'module_java']],
        name = "chromimpute_generate_train_data",
        command = """\
java {java_options} -jar /lb/project/mugqic/projects/rami_test/Tools/chromimpute/ChromImpute.jar \\
  GenerateTrainData \\
  -c {chrom} \\
  -r {resolution} \\
  {converteddir} \\
  {distancedir} \\
  {inputinfofile} \\
  {chrom_sizes} \\
  {output_dir} \\
  {mark}""".format(
        java_options = config.param('DEFAULT', 'java_options'),
        chrom = config.param('chromimpute','chrom'),
        resolution = config.param('chromimpute', 'resolution'),
        converteddir = converteddir,
        distancedir = input_dir,
        inputinfofile = config.param('chromimpute', 'inputinfofile'),
        chrom_sizes = config.param('chromimpute', 'chromsizes'),
        output_dir = output_dir,
        mark = mark
        )
    )


def train(input_dir, output_dir, sample, mark):
    return Job(
        [input_dir],
        [output_dir],
        [['java', 'module_java']],
        name = "chromimpute_train",
        command = """\
java {java_options} -jar /lb/project/mugqic/projects/rami_test/Tools/chromimpute/ChromImpute.jar \\
  Train \\
  {traindatadir} \\
  {inputinfofile} \\
  {predictordir} \\
  {sample} \\
  {mark}""".format(
        java_options = config.param('DEFAULT', 'java_options'),
        traindatadir = input_dir,
        inputinfofile = config.param('chromimpute', 'inputinfofile'),
        predictordir = output_dir,
        sample = sample,
        mark = mark
        )
    )

def apply(input_dir, converteddir, distancedir, predictordir, output_dir, sample, mark):
    return Job(
        [input_dir],
        [output_dir],
        [['java', 'module_java']],
        name = "chromimpute_apply",
        command = """\
java {java_options} -jar /lb/project/mugqic/projects/rami_test/Tools/chromimpute/ChromImpute.jar \\
    Apply \\
    -c {chrom} \\
    -r {resolution} \\
    {converteddir} \\
    {distancedir} \\
    {predictordir} \\
    {inputinfofile} \\
    {chrom_sizes} \\
    {output_dir} \\
    {sample} \\
    {mark}""".format(
        java_options = config.param('DEFAULT', 'java_options'),
        chrom = config.param('chromimpute','chrom'),
        resolution = config.param('chromimpute', 'resolution'),
        converteddir = converteddir,
        distancedir = distancedir,
        predictordir = predictordir,
        inputinfofile = config.param('chromimpute', 'inputinfofile'),
        chrom_sizes = config.param('chromimpute', 'chromsizes'),
        output_dir = output_dir,
        sample = sample,
        mark = mark
        )
    )

def eval(input_dir, percent1, percent2, converteddir, convertedFile, imputedFile, output_path):
    return Job(
        [input_dir],
        [],
        [['java', 'module_java']],
        name = "chromimpute_eval",
        command = """\
java {java_options} -jar /lb/project/mugqic/projects/rami_test/Tools/chromimpute/ChromImpute.jar \\
    Eval \\
    -p {percent1} {percent2} \\
    {converteddir} \\
    {convertedfile} \\
    {input_dir} \\
    {imputedfile} \\
    {chrom_sizes} > {output_path}""".format(
        java_options = config.param('DEFAULT', 'java_options'),
        percent1 = percent1,
        percent2 = percent2,
        converteddir = converteddir,
        convertedfile = convertedFile,
        input_dir = input_dir,
        imputedfile = imputedFile,
        chrom_sizes = config.param('chromimpute', 'chromsizes'),
        output_path = output_path
        )
    )



