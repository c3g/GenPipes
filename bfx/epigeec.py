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

def tohdf5(input_dir, output_dir, signal_file):
    signal_file_basename = os.path.basename(signal_file)+".hdf5"
    output_path = os.path.join(output_dir, signal_file_basename)

    return Job(
        [input_dir],
        [output_dir, output_dir+"_"+signal_file_basename], # Unique output name for each file
        [['epigeec', 'module_python']],
        name = "epigeec_tohdf5",
        command = """\
epigeec to_hdf5 \\
  -{file_type} \\
  {signal_file} \\
  {chromsizes} \\
  {resolution} \\
  {output}""".format(
            file_type = config.param('epigeec', 'file_type'),
            signal_file = signal_file,
            chromsizes = config.param('DEFAULT', 'chromsizes'),
            resolution = config.param('epigeec', 'resolution'),
            output = output_path
        )
    )

def filter(input_dir, output_dir, hdf5_file):
    hdf5_basename = os.path.basename(hdf5_file)

    return Job(
        [input_dir+"_"+hdf5_basename], # Job starts as soon as its corresponding bigwig file is converted to hdf5
        [output_dir],
        [['epigeec', 'module_python']],
        name = "epigeec_filter",
        command = """\
epigeec filter \\
  {hdf5_file} \\
  {chromsizes} \\
  {output_file} \\
  {select} \\
  {exclude}""".format(
        hdf5_file = hdf5_file,
        chromsizes = config.param('DEFAULT', 'chromsizes'),
        output_file = os.path.join(output_dir, hdf5_basename),
        select = config.param('epigeec', 'select'),
        exclude = config.param('epigeec', 'exclude')
        )
    )

def correlate(input_dir, output_dir, hdf5_list):
    return Job(
        [input_dir],
        [output_dir],
        [['epigeec', 'module_python']],
        name = "epigeec_correlate",
        command = """\
epigeec correlate \\
  {hdf5_list} \\
  {chromsizes} \\
  {outMatrix}""".format(
        hdf5_list = hdf5_list,
        chromsizes = config.param('DEFAULT', 'chromsizes'),
        outMatrix = os.path.join(output_dir, "correlation_matrix.tsv")
        )
    )

