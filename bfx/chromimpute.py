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

from core.config import *
from core.job import *

#chrom_sizes = os.path.join("../../ressources",config.param('chromimpute','chromsizes'))

def createInputInfo(path_to_dataset, output):
    filename = os.path.join(output, "inputinfo.txt")
    output = open(filename, "w+")
    for filename in os.listdir(path_to_dataset):
        parsed_file = filename.split(".")
        sample = parsed_file[2]
        mark = parsed_file[-3]
        output.write(sample+"\t"+mark+"\t"+filename+"\n")
    output.close()

def convert(self, path_to_dataset):
    createInputInfo(path_to_dataset, ".")
    #TEST command : REMOVE WHEN DONE
    command = """\
mkdir ChromImpute \\
java {java_options} ChromImpute.jar \\
  Convert \\
  -c {chrom} \\
  {path_to_dataset} \\
  inputinfo.txt \\
  {chrom_sizes} \\
  converteddir""".format(
        chrom = congif.param('chromimpute','chrom'),
        path_to_dataset = path_to_dataset,
        chrom_sizes = chrom_sizes
        )
    print(command)
    return Job(
        [path_to_dataset,"inputinfo.txt"],
        ['converteddir'],
        [
            ['chromimpute','module_java']
            ['chromimpute','module_chromimpute']
        ],
        command = """\
mkdir ChromImpute \\
java {java_options} ChromImpute.jar \\
  Convert \\
  -c {chrom} \\
  {path_to_dataset} \\
  inputinfo.txt \\
  {chrom_sizes} \\
  converteddir""".format(
        chrom = congif.param('chromimpute','chrom'),
        path_to_dataset = path_to_dataset,
        chrom_sizes = chrom_sizes
        )
    )

