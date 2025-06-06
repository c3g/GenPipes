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
import re

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def align(input1, input2, output_directory, outputs):

    inputs = []
    inputs.append(input1)
    if input2:
        inputs.append(input2)

    return Job(
        inputs,
        outputs,
        [
            ['bismark_align', 'module_bismark'],
            ['bismark_align', 'module_bowtie'],
            ['bismark_align', 'module_samtools']
        ],
        command="""\
bismark -q \\
  {other_options} \\
  {min_insert_size_PE} \\
  {genome_directory} \\
  {input1} \\
  {input2} \\
  --output_dir {output_directory} \\
  --temp_dir {tmp_dir}""".format(
            other_options=global_conf.global_get('bismark_align', 'other_options'),
            min_insert_size_PE=("-X " + global_conf.global_get('bismark_align', 'min_insert_size_PE')) if input2 else "",
            genome_directory=global_conf.global_get('bismark_align', 'bismark_assembly_dir'),
            input1="-1 "+input1 if input2 else input1,
            input2="-2 "+input2 if input2 else "",
            output_directory=output_directory,
            tmp_dir=global_conf.global_get('bismark_align', 'tmp_dir')
        )
    )


def dedup(input, outputs, library_type="PAIRED_END"):

    return Job(
        [input],
        outputs,
        [
            ['bismark_dedup', 'module_bismark'],
            ['bismark_dedup', 'module_bowtie'],
            ['bismark_dedup', 'module_samtools']

        ],
        command="""\
deduplicate_bismark \\
  {library} \\
  {other_options} \\
  {input}""".format(
            other_options=global_conf.global_get('bismark_dedup', 'other_options'),
            library="-p" if library_type=="PAIRED_END" else "-s",
            input=input
        ),
        removable_files=[re.sub(".bam", ".deduplicated.bam", input)]
    )

def methyl_call(input, outputs, library_type="PAIRED_END"):

    return Job(
        [input],
        outputs,
        [
            ['bismark_methyl_call', 'module_bismark'],
            ['bismark_methyl_call', 'module_samtools']

        ],
        command="""\
bismark_methylation_extractor \\
  {library} \\
  {other_options} \\
  --output {output_directory} \\
  {input}""".format(
            other_options=global_conf.global_get('bismark_methyl_call', 'other_options'),
            library="-p" if library_type=="PAIRED_END" else "-s",
            input=input,
            output_directory=os.path.dirname(outputs[0])
        )
    )

def bed_graph(inputs, output_prefixe, output_directory):

    if not isinstance(inputs, list):
        inputs=[inputs]
    outputs = [
        os.path.join(output_directory, output_prefixe + ".bedGraph.gz"),
        os.path.join(output_directory, output_prefixe + ".bismark.cov.gz")
    ]

    return Job(
        inputs,
        outputs,
        [
            ['bismark_bed_graph', 'module_bismark']
        ],
        command="""\
bismark2bedGraph \\
  {other_options} \\
  --output {output} \\
  --dir {directory} \\
  {inputs}""".format(
            inputs="".join([" \\\n  " + input for input in inputs]),
            output=output_prefixe + ".bedGraph.gz",
            directory=output_directory,
            other_options=global_conf.global_get('bismark_bed_graph', 'other_options')
        )
    )

def coverage2cytosine(input, output):
    return Job(
        [input],
        [output],
        [
            ['bismark_coverage2cytosine', 'module_bismark']
        ],
        command="""\
coverage2cytosine \\
  {other_options} \\
  --output {output} \\
  {input}""".format(
            input=input,
            output=output,
            other_options=global_conf.global_get('bismark_coverage2cytosine', 'other_options')
        )
    )
