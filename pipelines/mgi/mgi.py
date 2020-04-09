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
import logging
import math
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs
from pipelines import common

from bfx import adapters
from bfx import cutadapt
from pipelines.dnaseq import dnaseq

log = logging.getLogger(__name__)


    

class MGISeq(dnaseq.DnaSeq):
    """
    MGI-Seq Pipeline
    ================

    pwet
    """

    def __init__(self, protocol=None):
        self._protocol=protocol
        # Add pipeline specific arguments
        super(MGISeq, self).__init__(protocol)

    
    def cutadapt(self):
        """
        """

        jobs = []

        for readset in self.readsets:
            output_dir = os.path.join("trim", readset.sample.name)
            trim_file_prefix = os.path.join(output_dir, readset.name)

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[readset.fastq1, readset.fastq2]]
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                adapter_fwd = readset.adapter1
                adapter_rev = readset.adapter2

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[readset.fastq1]]
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                adapter_fwd = readset.adapter1
                adapter_rev = None

            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            jobs.append(
                concat_jobs([
                    Job(
                        command="mkdir -p " + output_dir,
                        removable_files=[output_dir],
                        samples=[readset.sample]
                        ),
                    cutadapt.trim(
                        fastq1,
                        fastq2,
                        trim_file_prefix,
                        adapter_fwd,
                        adapter_rev
                        )
                    ], name="cutadapt_trimming." + readset.name)
                )

        return jobs

    

    @property
    def steps(self):
        return [
            self.cutadapt,
            self.bwa_mem_sambamba_sort_sam,
            self.sambamba_merge_sam_files,
        ]

if __name__ == '__main__':
    MGISeq()
