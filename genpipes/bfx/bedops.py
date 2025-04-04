################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
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

# Python Standard Modules
import os

# MUGQIC Modules
from ..core.job import Job


def sort_bed(input_bed, output_file, other_options=""):
    return Job(
        [input_bed],
        [output_file],
        [
            ['bed_sort', 'module_bedops']
        ],
        name="bed_sort." + os.path.basename(input_bed),
        command="""
sort-bed {other_options} \\
  {input} \\
  > {output}""".format(
            other_options=other_options,
            input=input_bed,
            output=output_file
        )
    )

def bedmap_echoMapId(bed1, bed2, output, split=True, other_options=""):
    return Job(
        [bed1, bed2],
        [output],
        [
            ['bedmap_echoMapId', 'module_bedops']
        ],
        name="bedmap_echoMapId." + os.path.basename(output),
        command="""
bedmap {other_options} \\
  --echo --echo-map-id \\
  {bed1} \\
  {bed2} \\
  {split} > {output}""" .format(
            other_options=other_options,
            bed1=bed1,
            bed2=bed2,
            split=" | tr '|' '\\t' " if split else "",
            output=output
        )
    )
