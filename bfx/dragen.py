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
                      input_dependency=None, output_dependency=None, protocol="hybrid"):
   # print(protocol)
    if protocol == "dragen":

        duplicate_marking = config.param('dragen_align', 'duplicate_marking', param_type='boolean')
        remove_duplicates = config.param('dragen_align', 'remove_duplicates', param_type='boolean')
    else:
        duplicate_marking = "false"
        remove_duplicates = "false"

    if input_dependency is not None:
        inputs = input_dependency
    else:
        if fastq2 is not None:
            inputs = [fastq1, fastq2]
        else:
            inputs = [fastq1]
    if output_dependency is not None:
        outputs = output_dependency
   # if duplicate_marking == "true":
        ##outputs = [os.path.join(output_dir, readsetName + ".sorted.bam")]
        output_prefix = readsetName + ".sorted"
    #else:
     #   output_prefix = readsetName
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
    --output-directory {output_dir} \\
    --output-file-prefix {output_prefix} \\
    --methylation-generate-cytosine-report {ct_report} \\
    --methylation-mapping-implementation {mapping_implementation} \\
    --enable-sort {sort} \\
    --enable-duplicate-marking {duplicate_marking} \\
    --RGID {rgid} \\
    --RGLB {rglb} \\
    --RGSM {rgsm} \\
    -1 {fastq1} \\
    {fastq2} \\
    --remove-duplicates {remove_duplicates} {other_options}""".format(
            output_dir=output_dir,
            output_prefix=output_prefix,
            dragen_tmp=config.param('dragen_align', 'tmp_dragen', param_type='string'),
            met_protocol=config.param('dragen_align', 'methylation_protocol', param_type='string'),
            reference=config.param('dragen_align', 'reference', param_type='string'),
            ct_report="true" if config.param('dragen_align', 'CTreport', param_type='boolean') else "false",
           # sort=config.param('dragen_align', 'sort', param_type='string'),
            sort="true",
            mapping_implementation=config.param('dragen_align', 'mapping_implementation', param_type='string'),
            duplicate_marking=duplicate_marking,
            remove_duplicates=remove_duplicates,
            rgid=readGroupID,
            rglb=libraryName,
            rgsm=sampleName,
            other_options=config.param('dragen_align', 'other_options', param_type='string', required=False),
            fastq1=fastq1,
            fastq2="-2 " + fastq2 if fastq2 != None else ""
        ),
    )


def call_methylation(bam, output_dir, readsetName,
                      input_dependency=None, output_dependency=None, protocol="dragen"):

    if input_dependency is not None:
        inputs = input_dependency
    else:
        inputs = [bam]
    if output_dependency is not None:
        outputs = output_dependency
        output_prefix = readsetName + ".sorted"
    else:
        #TODO complete this to get output file path
        outputs = []
    return Job(
        inputs,
        outputs,
        [],
        command="""\
dragen_reset && \\
dragen --enable-methylation-calling true \\
    --intermediate-results-dir {dragen_tmp} \\
    --ref-dir {reference} \\
    --output-directory {output_dir} \\
    --output-file-prefix {output_prefix} \\
    --methylation-reports-only {ct_report} \\
    --methylation-generate-cytosine-report {ct_report} \\
    --methylation-generate-mbias-report {mbias_report} \\
    --enable-sort {sort} \\
    --methylation-match-bismark {match_bismark} \\
    -b {bam} {other_options}""".format(
            output_dir=output_dir,
            output_prefix=output_prefix,
            dragen_tmp=config.param('dragen_align', 'tmp_dragen', param_type='string'),
            reference=config.param('dragen_align', 'reference', param_type='string'),
            ct_report="true",
            sort="false",
            mbias_report="true",
            match_bismark="true",
            mapping_implementation=config.param('dragen_align', 'mapping_implementation', param_type='string'),
            other_options=config.param('dragen_align', 'other_options', param_type='string', required=False),
            bam=bam
        ),
    )

