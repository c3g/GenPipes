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

# !/usr/bin/env python

# Python Standard Modules
import logging
import math
import os
import re
import sys
import itertools

# MUGQIC Modules
# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs
from bfx.readset import parse_illumina_readset_file


def align_methylation(fastq1, fastq2, output_dir, readsetName, sampleName, libraryName, readGroupID,
                      input_dependency=None, output_dependency=None):
    if input_dependency is not None:
        inputs = input_dependency
    else:
        if fastq2 is not None:
            inputs = [fastq1, fastq2]
        else:
            inputs = [fastq1]
    if output_dependency is not None:
        outputs = output_dependency
    else:
        outputs = [os.path.join(output_dir, readsetName + ".bam")]

    return Job(
        inputs,
        outputs,
        [],
        command="""\
dragen_reset && \\
dragen --enable-methylation-calling true \\
    --intermediate-results-dir {dragen_tmp} \\
    --methylation-protocol {met_protocol} \\
    --ref-dir {reference} \\
    --output-directory  {output_dir} \\
    --output-file-prefix {readset} \\
    --methylation-generate-cytosine-report {ct_report} \\
    --enable-sort false \\
    --RGID {rgid} \\
    --RGLB {rglb} \\
    --RGSM {rgsm} \\
    -1 {fastq1} \\
    {fastq2} {other_options}""".format(
            output_dir=output_dir,
            readset=readsetName,
            dragen_tmp=config.param('dragen_align', 'tmp_dragen', type='string'),
            met_protocol=config.param('dragen_align', 'methylation_protocol', type='string'),
            reference=config.param('dragen_align', 'reference', type='string'),
            ct_report="true" if config.param('dragen_align', 'CTreport', type='boolean') else "false",
            rgid=readGroupID,
            rglb=libraryName,
            rgsm=sampleName,
            other_options=config.param('dragen_align', 'other_options', type='string'),
            fastq1=fastq1,
            fastq2="-2 " + fastq2 if fastq2 != None else ""
        ),
    )

