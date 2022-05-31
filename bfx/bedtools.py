################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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
from core.config import *
from core.job import *

from bfx import ucsc

def graph(input_bam, output_bed_graph, library_type="PAIRED_END"):

    if library_type == "PAIRED_END":
        if "forward" in output_bed_graph:
            samtools_options = "-F 256 -f 81 "
        elif "reverse" in output_bed_graph:
            samtools_options = "-F 256 -f 97 "
        else:
            raise Exception("Error: PAIRED_END library was provided but no strand orientation could be determined from " + output_bed_graph + "...")
    else:
        samtools_options = "-F 256"

    return Job(
        [input_bam],
        [output_bed_graph],
        [
            ['bedtools', 'module_samtools'],
            ['bedtools', 'module_bedtools']
        ],
        command="""\
nmblines=$(samtools view {samtools_options} {input_bam} | wc -l) && \\
if [[ $nmblines -ge 100000 ]]
  then
    scalefactor=$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc)
  else
    scalefactor=1
fi
genomeCoverageBed {other_options} -bg -split -scale $scalefactor \\
  -ibam {input_bam} \\
  -g {chromosome_size} \\
  > {output_bed_graph}""".format(
            samtools_options=samtools_options,
            input_bam=input_bam,
            chromosome_size=config.param('bedtools_graph', 'chromosome_size', param_type='filepath'),
            other_options=config.param('bedtools_graph', 'other_options', required=False),
            output_bed_graph=output_bed_graph
        )
    )

def intersect(input_bam, output_bam, target_bed, include_header=False):

    return Job(
        [input_bam],
        [output_bam],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools intersect {other_options} {include_header} \\
  -a {input_bam} \\
  -b {target_bed} \\
  > {output_bam}""".format(
            input_bam=input_bam,
            target_bed=target_bed,
            other_options=config.param('bedtools_intersect', 'other_options', required=False),
            include_header="-header" if include_header else "",
            output_bam=output_bam
        )
    )

def intersect_beds(bed1, bed2, output_bed, other_options=""):

    return Job(
        [bed1, bed2],
        [output_bed],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools intersect {other_options} \\
  -a {bed1} \\
  -b {bed2} \\
  > {output_bed}""".format(
            bed1 = bed1,
            bed2 = bed2,
            other_options = other_options,
            output_bed = output_bed
        )
    )

def bamtobed(input_bam, output_bed):

    return Job(
        [input_bam],
        [output_bed],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools bamtobed {other_options} \\
  -i {input_bam}{output_bed}""".format(
            input_bam=input_bam,
            other_options=config.param('bedtools_coverage', 'other_options', required=False),
            output_bed=" \\\n  > " + output_bed if output_bed else ""
        )
    )

def coverage(input_file, output_file):

    return Job(
        [input_file],
        [output_file],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools coverage {other_options} \\
  -a {intervals} \\
  -b {input} \\
  > {output_file}""".format(
            intervals=config.param('bedtools_coverage', 'gc_intervals'),
            input=input_file,
            other_options=config.param('bedtools_coverage', 'other_options', required=False),
            output_file=output_file
        )
    )


def genomecov(input_file, output_file):

    return Job(
        [input_file],
        [output_file],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools genomecov {other_options} \\
  {input} \\
  {genome} \\
  > {output_file}""".format(
            input="-ibam " + input_file if re.search("\.bam$", os.path.basename(input_file)) else "-i " + input_file,
            genome="-g " + config.param('bedtools_genomecov', 'genome_fasta') if not re.search("\.bam$", os.path.basename(input_file)) else "",
            other_options=config.param('bedtools_genomecov', 'other_options', required=False),
            output_file=output_file
        )
    )

def bamtofastq(input_bam, output_pair1, output_pair2, other_options=config.param('bedtools_bamtofastq', 'other_options', required=False), pigz_threads=config.param('bedtools_bamtofastq', 'pigz_threads', required=False)):
    if output_pair2:  # Paired end reads
        outputs = [output_pair1, output_pair2]
    else:   # Single end reads
        outputs = [output_pair1]

    return Job(
        [input_bam],
        [outf + ".gz" for outf in outputs],
        [
            ['bedtools', 'module_bedtools'],
            ['pigz', 'module_pigz']
        ],
        command="""\
bedtools bamtofastq {other_options} \\
  -i {input_bam} \\
  {output_pair1} \\
  {output_pair2} && \\
pigz -f -p {pigz_threads} {input_fq}""".format(
    input_bam=input_bam,
    other_options=other_options,
    output_pair1="-fq " + output_pair1,
    output_pair2="-fq2 " + output_pair2 if output_pair2 else "",
    pigz_threads=pigz_threads,
    input_fq=" ".join(outputs)
    )
        )
