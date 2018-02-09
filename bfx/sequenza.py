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


from core.config import *
from core.job import *

#Assuming that the module file has already been written and the scripts have been added to path

def genome_gz(output):
    return Job(
        [None],
        [output],
        [
            ['sequenza', 'module_python'],
            ['sequenza', 'module_sequenza'],
        ],
        command="""\\
python $SEQUENZA_BIN/sequenza-utils.py \\
    GC-windows -w {window} \\
    {input} \\
    | gzip > \\
    {out}""".format(
        input=config.param('samtools_mpileup','genome_fasta',type='filepath'),
        out=output,
        window=config.param('sequenza','window_length')
    ),  
  )


def sequenza_seqz(normal_gz, tumor_gz, genome, output):
    return Job(
        [normal_gz, tumor_gz],
        [output],
        [
            ['sequenza', 'module_python'],
            ['sequenza', 'module_sequenza'],
        ],
        command="""\\
python $SEQUENZA_BIN/sequenza-utils.py \\
    pileup2seqz  {pileup_options} \\
    -gc  {gen}   \\
    -n {normal}  \\
    -t {tumor}   \\
    {out}""".format(
        gen=genome,
        normal=normal_gz,
        tumor=tumor_gz,
        pileup_options=config.param('sequenza','pileup_options'),
        out=" \\\n > " + output if output else ""
        )
    )

def sequenza_bin(seqz_gz, output):
    return Job(
        [seqz_gz],
        [output],
        [
            ['sequenza', 'module_python'],
            ['sequenza', 'module_sequenza'],
        ],
        command="""\\
python $SEQUENZA_BIN/sequenza-utils.py  \\
    seqz-binning  \\
    -w {window}  \\
    -s {seqz_gz} \\
    {output}""".format(
        window=config.param('sequenza','bin_window_size'),
        seqz_gz=seqz_gz,
        output=" \\\n > " + output if output else "",
        )
    )

def sequenza_main(seqz,output_folder,sample_name):
    return Job(
        [seqz],
        [None],
        [
            ['sequenza', 'module_mugqic_tools'],
            ['sequenza', 'module_R'],
        ],
        command="""\\
Rscript $R_TOOLS/RunSequenza_analysis.R \\
    {input}   \\
    {OUTPUT_FOLDER}   \\
    {OUTPUT_BASE_NAME}""".format(
            input=seqz,
            OUTPUT_FOLDER=output_folder,
            OUTPUT_BASE_NAME=sample_name
        )
    )
