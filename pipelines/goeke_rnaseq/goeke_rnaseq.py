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
import argparse
import logging
import math
import os
import re
import sys
import subprocess

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs
import utils.utils

from bfx import adapters
from bfx import fastqc
from bfx import salmon
from bfx import bash_cmd as bash

from pipelines import common
from pipelines.rnaseq import rnaseq

log = logging.getLogger(__name__)

class RnaSeqLight(rnaseq.RnaSeq):
    def __init__(self,protocol=None):
        self._protocol=protocol
        super(RnaSeqLight, self).__init__(protocol)

    def fastqc(self):
        """
        Step 1: Quality Control (with FastQC)
        """
        jobs = []

        for readset in self.readsets:

            output_dir = os.path.join(self.output_dir, "qc", readset.sample.name, readset.name)
            job_name = "fastqc." + readset.name

            adapter_file = config.param('fastqc', 'adapter_file', required=True, type='filepath')

            # PAIRED
            if readset.run_type == "PAIRED_END":
                fastq1 = readset.fastq1
                fastq2 = readset.fastq2
                output = os.path.join(output_dir, readset.name + "_fastqc.zip")
                job = fastqc.fastqc(input1=fastq1, input2=fastq2, output_directory=output_dir,
                                    output=output, adapter_file=adapter_file)
                job_samples = [readset.sample]

                jobs.append(
                    concat_jobs([bash.mkdir(output_dir, remove=True), job], name=job_name, samples=job_samples)
                )

            # SINGLE
            elif readset.run_type == "SINGLE_END":
                fastq1 = readset.fastq1
                fastq2 = None
                output = os.path.join(output_dir, readset.name + "_fastqc.zip")
                job = fastqc.fastqc(input1=fastq1, input2=fastq2, output_directory=output_dir,
                                    output=output, adapter_file=adapter_file)
                job_samples = [readset.sample]

                jobs.append(
                    concat_jobs([bash.mkdir(output_dir, remove=True), job], name=job_name, samples=job_samples)
                )

        return jobs


    def salmon_index(self):
        """
        Step 2: Index Creation (with Salmon)
        """
        jobs = []
        output_dir = os.path.join(self.output_dir, "salmon_index")
        transcriptome_file = config.param('salmon_index', 'transcriptome_fasta', required=True, type='filepath')

        jobs.append(concat_jobs([bash.mkdir(output_dir, remove=True),
                                 salmon.salmon_index(transcriptome_file, output_dir)],
                                 name="salmon_index"))

        return jobs

    def salmon_quant(self):
        """
        Step 3: Quantification (with Salmon)
        """
        jobs = []

        for readset in self.readsets:

            output_dir = os.path.join(self.output_dir, "salmon", readset.sample.name, readset.name)
            salmon_idx= os.path.join(self.output_dir, "salmon_index")
            job_name = "salmon_quant." + readset.name

            # PAIRED
            if readset.run_type == "PAIRED_END":
                fastq1 = readset.fastq1
                fastq2 = readset.fastq2
                job = salmon.salmon_quant(readset_name=readset.name, read1_fastq=fastq1, read2_fastq=fastq2,
                                          output_directory=output_dir, salmon_index=salmon_idx)
                job_samples = [readset.sample]

                jobs.append(
                    concat_jobs([bash.mkdir(output_dir, remove=True), job], name=job_name, samples=job_samples)
                )

            # SINGLE
            elif readset.run_type == "SINGLE_END":
                fastq1 = readset.fastq1
                fastq2 = None
                job = salmon.salmon_quant(readset_name=readset.name, read1_fastq=fastq1, read2_fastq=fastq2,
                                          output_directory=output_dir, salmon_index=salmon_idx)
                job_samples = [readset.sample]

                jobs.append(
                    concat_jobs([bash.mkdir(output_dir, remove=True), job], name=job_name, samples=job_samples)
                )

        return jobs


############

    @property
    def steps(self):
        return [
            self.fastqc,
            self.salmon_index,
            self.salmon_quant,
            ]

if __name__ == '__main__':
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        RnaSeqLight()

