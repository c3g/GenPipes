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
import commands
import gzip
import subprocess
import pysam

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.design import *
from pipelines import common
from pipelines.hicseq import hicseq


from bfx import picard
from bfx import samtools
from bfx import hicup
from bfx import homer
from bfx import multiqc
from bfx import genome
from bfx import bedtools
from bfx import chicago

log = logging.getLogger(__name__)

class ChicSeq(hicseq.HicSeq):
    """
    Capture Hi-C Pipeline
    =====================

    Hi-C experiments allow researchers to understand chromosomal folding and structure using proximity ligation techniques.
    Capture Hi-C experiments allow researchers to focus on regions of interest by capturing those regions on a chip.
    The pipeline starts by trimming adaptors and low quality bases. It then maps the reads to a reference genome using HiCUP.
    HiCUP first truncates chimeric reads, maps reads to the genome, then it filters Hi-C artifacts and removes read duplicates.
    Samples from different lanes are merged. CHiCAGO is then used to filter capture-specific artifacts and call significant 
    interactions. This pipeline also identifies enrichement of regulatory features when provided with ChIP-Seq marks and also 
    identifies enrichement of GWAS SNPs.
    """

    def __init__(self):
        super(ChicSeq, self).__init__()



    @property
    def output_dirs(self):
        dirs = {'hicup_output_directory': 'alignment',
                'bams_output_directory': 'alignment',
                'chicago_input_files': 'input_files',
                'chicago_output_directory': 'chicago'
                }

        return dirs

    def create_rmap_file(self):
        ## return 1 rmap per enzyme
        output = os.path.join(self.output_dirs['chicago_input_files'], self.enzyme + ".rmap")
        input_file = self.genome_digest

        command = """mkdir -p {dir} && \\
        cut -f 1-3 {input_file} > {output}.tmp && \\
        awk 'BEGIN {{FS=\"\\t\"; OFS=\"\\t\"}} NR>2 {{print $0, NR}}' {output}.tmp > {output} && \\
        rm {output}.tmp""".format(dir = self.output_dirs['chicago_input_files'], input_file = input_file, output = output)

        return [Job(input_files = [input_file], 
            output_files = [output], 
            command = command, 
            name = "create_rmap_file." + self.enzyme)]


    def create_baitmap_file(self):
        ## return 1 baitmap per enzyme/capture array combination
        
        input_rmap = os.path.join(self.output_dirs['chicago_input_files'], self.enzyme + ".rmap")
        input_bait = config.param('create_baitmap_file', "baitBed")
        output_file = re.sub("\.bed", "", os.path.basename(input_bait)) + "_" + self.enzyme + ".baitmap"
        output = os.path.join(self.output_dirs['chicago_input_files'], output_file)
        annotation = config.param('create_baitmap_file', "annotation")

        job_intersectBeds = bedtools.intersect_beds(input_rmap, input_bait, output + ".tmp", "-wa -u")

        job_anno = Job(input_files = [output + ".tmp"], 
            output_files = [output], 
            command = """awk 'BEGIN {{FS=\"\\t\"; OFS=\"\\t\"}}{{print $0, \"{annotation}\"NR}}' {outputTmp}  >  {output}""".format(annotation = annotation, outputTmp = output + ".tmp", output = output),
            name = "create_baitmap_file.addAnno." + output,
            removable_files = [output + ".tmp"]
            )

        job = concat_jobs([job_intersectBeds, job_anno])
        job.name = "create_baitmap_file." + output_file
        return [job]


    def create_design_files(self):
        rmapfile = os.path.join(self.output_dirs['chicago_input_files'], self.enzyme + ".rmap")
        baitmapfile = os.path.join(self.output_dirs['chicago_input_files'], os.path.basename(re.sub("\.bed$", "", config.param('create_baitmap_file', "baitBed")) + "_" + self.enzyme + ".baitmap"))
        other_options = config.param('create_design_files', 'other_options')
        designDir = self.output_dirs['chicago_input_files']
        outfilePrefix = os.path.join(self.output_dirs['chicago_input_files'], os.path.basename(re.sub("\.bed$", "", config.param('create_baitmap_file', "baitBed")) + "_" + self.enzyme))
        job = chicago.makeDesignFiles(rmapfile, baitmapfile, outfilePrefix, designDir, other_options)

        return [job]


    def create_input_files(self):
        jobs = []
        rmapfile = os.path.join(self.output_dirs['chicago_input_files'], self.enzyme + ".rmap")
        baitmapfile = os.path.join(self.output_dirs['chicago_input_files'], os.path.basename(re.sub("\.bed$", "", config.param('create_baitmap_file', "baitBed")) + "_" + self.enzyme + ".baitmap"))
        other_options=""


        for sample in self.samples:
            name = os.path.join(self.output_dirs['chicago_input_files'], sample.name)
            bam = os.path.join(self.output_dirs['bams_output_directory'], sample.name, sample.name + ".merged.bam")
            job = chicago.bam2chicago(bam, baitmapfile, rmapfile, name, other_options)
            jobs.append(job)

        return jobs


    @property
    def steps(self):
        return [
            self.samtools_bam_sort,
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.fastq_readName_Edit,
            self.hicup_align,
            self.samtools_merge_bams,
            self.create_rmap_file,
            self.create_baitmap_file,
            self.create_design_files,
            self.create_input_files,
            self.multiqc_report
        ]

if __name__ == '__main__':
    ChicSeq()
