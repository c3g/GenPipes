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

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def verify(input_bam, output_prefix, bed_file=None):
    if not isinstance(input_bam, list):
        inputs = [input_bam]
        
    if bed_file is not None:
        inputs.append(bed_file)
        
    return Job(
        inputs,
        [output_prefix + ".selfSM"],
        [
            ['verify_bam_id2', 'module_verify_bam_id2']
        ],
        command="""\
VerifyBamID {other_options} \\
  --SVDPrefix {svdprefix} \\
  --Reference {reference} \\
  --BamFile {input_bam} \\
  --Output {output_prefix}{bed_file}""".format(
            svdprefix=config.param('verify_bam_id2', 'svd_dataset'),
            reference=config.param('verify_bam_id2', 'genome_fasta', param_type='filepath'),
            input_bam=input_bam,
            output_prefix=output_prefix,
            bed_file=" --BedPath " + bed_file if bed_file else "",
            other_options=config.param('verify_bam_id2', 'options')
        )
    )

def parse_contamination_freemix_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export contamination_freemix=`grep -v "^#" {input_file} | cut -f7`"""
        )