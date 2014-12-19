#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def graph(input_bam, output_bed_graph, output_wiggle):

    return Job(
        [input_bam],
        [output_bed_graph, output_wiggle],
        [
            ['bedtools', 'module_samtools'],
            ['bedtools', 'module_bedtools'],
            ['bedtools', 'module_ucsc']
        ],
        command="""\
nmblines=$(samtools view -F 256 -f 81 {input_bam} | wc -l) && \\
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \\
genomeCoverageBed -bg -split -scale $scalefactor \\
  -ibam {input_bam} \\
  -g {chromosome_size} \\
  > {output_bed_graph} && \\
bedGraphToBigWig \\
  {output_bed_graph} \\
  {chromosome_size} \\
  {output_wiggle}""".format(
        input_bam=input_bam,
        chromosome_size=config.param('bedtools', 'chromosome_size', type='filepath'),
        output_bed_graph=output_bed_graph,
        output_wiggle=output_wiggle
        )
    )
