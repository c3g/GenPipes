################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import os

from ..core.config import global_conf
from ..core.job import Job

def tohdf5(output_dir, signal_file, ini_section='epigeec'):
    signal_file_basename = os.path.basename(signal_file)+".hdf5"
    output_path = os.path.join(output_dir, signal_file_basename)
    log.info(signal_file)
    return Job(
        [signal_file],
        [output_path], # Unique output name for each file
        [
            [ini_section, 'module_python']
        ],
        command = """\
epigeec to_hdf5 \\
  -{file_type} \\
  {signal_file} \\
  {chromsizes} \\
  {resolution} \\
  {output}""".format(
            file_type = global_conf.global_get(ini_section, 'file_type'),
            signal_file = signal_file,
            chromsizes = global_conf.global_get(ini_section, 'chromosome_size'),
            resolution = global_conf.global_get(ini_section, 'resolution'),
            output = output_path
        )
    )

def filter(output_dir, hdf5_file, ini_section='epigeec'):

##linux if command below checks whether select and exclude files are available. then it adds -s and -e options to the
# command based on their avaialability
    filtered_output = os.path.join(output_dir, os.path.basename(hdf5_file))
    return Job(
        [hdf5_file], # Job starts as soon as its corresponding bigwig file is converted to hdf5
        [filtered_output],
        [
            [ini_section, 'module_python']
        ],
        command = """\
epigeec filter \\
  {hdf5_file} \\
  {chromsizes} \\
  {output_file} \\
  `if [[ "{select}" != "" ]]; then echo "-s {select}"; fi; if [[ "{exclude}" != "" ]]; then echo "-e {exclude}" ;fi`""".format(
        hdf5_file = hdf5_file,
        chromsizes = global_conf.global_get(ini_section, 'chromosome_size'),
        output_file = filtered_output,
        select = global_conf.global_get(ini_section, 'select'),
        exclude = global_conf.global_get(ini_section, 'exclude')
        )
    )

def correlate(input_files, output_dir, hdf5_list, ini_section='epigeec'):
    correlation_matrix = os.path.join(output_dir, "correlation_matrix.tsv")
    return Job(
        input_files,
        [correlation_matrix],
        [
            [ini_section, 'module_python']
        ],
        name = "epigeec_correlate",
        command = """\
epigeec correlate \\
  {hdf5_list} \\
  {chromsizes} \\
  {outMatrix}""".format(
        hdf5_list = hdf5_list,
        chromsizes = global_conf.global_get(ini_section, 'chromosome_size'),
        outMatrix = correlation_matrix
        )
    )

def generate_hdf5_list(hdf5_list_file, file_path):
    return Job(
        input_files=[hdf5_list_file],
        command="""\
            echo -e "{file_path}" >> {hdf5_list_file}""".format(
            hdf5_list_file=hdf5_list_file,
            file_path=file_path
            )
        )

