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
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.design import *
from pipelines import common

from bfx import picard
from bfx import samtools

log = logging.getLogger(__name__)

class HicSeq(common.Illumina):
    """
    Hi-C Pipeline
    ==============

    Hi-C experiments allow researchers to understand chromosomal folding and structure using proximity ligation techniques.
    The pipeline starts by trimming adaptors and low quality bases. It then maps the reads to a reference genome using HiCUP.
    HiCUP first truncates chimeric reads, maps reads to the genome, then it filters Hi-C artifacts and removes read duplicates.
    Samples from different lanes are merged and a tag directory is created by Homer, which is also used to produce the interaction
    matrices and compartments. TopDom is used to predict topologically associating domains (TADs) and homer is used to identify
    significant interactions.

    An example of the Hi-C report for an analysis on public data (GM12878 Rao. et al.) is available for illustration purpose only:
    [Hi-C report](<url>).

    [Here](<url>) is more information about Hi-C pipeline that you may find interesting.
    """

    def __init__(self):
        super(HicSeq, self).__init__()

    #@property
    def hicup_align(self):
        """
        Paired-end Hi-C reads are truncated, mapped and filtered using HiCUP. The resulting bam file is filtered for Hi-C artifacts and
        duplicated reads. It is ready for use as input for downstream analysis.

        For more detailed information about the HICUP process visit: [HiCUP] (https://www.bioinformatics.babraham.ac.uk/projects/hicup/overview/)
        """


        jobs = []

        # create hicup output directory:
        output_directory = "HiCUP_Alignments"


        for readset in self.readsets:
            sample_output_dir = os.path.join(output_directory, readset.name)
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")

            if readset.run_type != "PAIRED_END":
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END for Hi-C analysis)!")

            candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
            if readset.fastq1 and readset.fastq2:
                candidate_input_files.append([readset.fastq1, readset.fastq2])
            if readset.bam:
                candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
            [fastq1, fastq2] = self.select_input_files(candidate_input_files)


            ## create HiCUP configuration file:
            configFileContent = """
            Outdir: {sample_output_dir}
            Threads: {threads}
            Quiet:{Quiet}
            Keep:{Keep}
            Zip:{Zip}
            Bowtie2: {Bowtie2_path}
            R: {R_path}
            Index: {Genome_Index_hicup}
            Digest: {Genome_Digest}
            Format: {Format}
            Longest: {Longest}
            Shortest: {Shortest}
            {fastq1}
            {fastq2}
            """.format(outDir = sample_output_dir,
                threads = config.param('hicup_align', 'threads'),
                Quiet = config.param('hicup_align', 'Quiet'),
                Keep = config.param('hicup_align', 'Keep'),
                Zip = config.param('hicup_align', 'Zip'),
                Bowtie2_path = config.param('hicup_align', 'Bowtie2_path'),
                R_path = config.param('hicup_align', 'R_path'),
                Genome_Index_hicup = config.param('hicup_align', 'Genome_Index_hicup'),
                Genome_Digest = config.param('hicup_align', 'Genome_Digest'),
                Format = config.param('hicup_align', 'Format'),
                Longest = config.param('hicup_align', 'Longest'),
                Shortest = config.param('hicup_align', 'Shortest'),
                fastq1 = fastq1,
                fastq2 = fastq2)

            ## write configFileContent to temporary file:
            with open("hicup_align." + readset.name + ".conf", "w") as conf_file:
                conf_file.write(configFileContent)

            ## hicup command
            command="mkdir -p {sample_output_dir} && hicup -c {conf_file}".format(sample_output_dir= sample_output_dir, conf_file=conf_file)
            job = Job(input_files= [fastq1, fastq2, conf_file],
                    output_files=[".".join(fastq1, fastq2, "hicup.bam")],
                    module_entries= ["module_bowtie2", "module_mugqic_R_packages", "module_perl", "module_R"],
                    name= "hicup_align." + readset.name,
                    command=command,
                    removable_files=[conf_file]
                    )

            jobs.append(job)
            return jobs



    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.hicup_align
        ]

if __name__ == '__main__':
    HicSeq()
