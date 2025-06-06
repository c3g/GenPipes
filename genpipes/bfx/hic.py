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
from ..core.config import global_conf
from ..core.job import Job

def create_input(bam_input, sample_name):
    return Job(
        [bam_input],
        [sample_name + ".juicebox.input", sample_name + ".juicebox.input.sorted"],
        [
            ["create_hic_file", "module_mugqic_tools"]
        ],
        command="""
bash CreateHicFileInput.sh \\
    {bam} \\
    {name} \\
    {tmpDir}""".format(
            bam=bam_input,
            name=sample_name,
            tmpDir=os.path.expandvars("$(pwd)")
        ),
        removable_files=[sample_name + ".juicebox.input", sample_name + ".juicebox.input.sorted", bam_input]
    )

def create_hic(juicebox_input, hic_output, assembly):
    return Job(
        [juicebox_input],
        [hic_output],
        [
            ["create_hic_file", "module_java"],
            ["create_hic_file", "module_juicer"]
        ],
        command="""
java -jar $juicer_JAR \\
  pre \\
  -q {q} \\
  {input} \\
  {output} \\
  {assembly}""".format(
	        q=global_conf.global_get('create_hic_file', 'q'),
            input=juicebox_input,
            output=hic_output,
            assembly=assembly
        )
    )
