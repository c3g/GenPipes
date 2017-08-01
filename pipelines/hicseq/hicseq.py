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
        self.argparser.add_argument("-e", "--enzyme", help = "Restriction Enzyme used to generate Hi-C library", choices = ["DpnII", "HindIII", "NcoI", "MboI"])
        super(HicSeq, self).__init__()


    @property
    def enzyme(self):
        return self.args.enzyme



    @property
    def restriction_site(self):
        """ sets the restriction enzyme recogntition site and genome digest location based on enzyme"""
        if (self.enzyme == "DpnII") or (self.enzyme == "MboI"):
            restriction_site = "GATC"
        elif self.enzyme == "HindIII":
            restriction_site = "AAGCTT"
        elif self.enzyme == "NcoI":
            restriction_site = "CCATGG"
        else:
            raise Exception("Error: Selected Enzyme is not yet available for Hi-C analysis!")
        return restriction_site

    @property
    def genome_digest(self):
        genome_digest = os.path.expandvars(config.param('hicup_align', "genome_digest_" + self.enzyme))
        return genome_digest

    def fastq_readName_Edit(self):
        """
        Removes the added /1 and /2 by picard's sam_to_fastq transformation to avoid issues with downstream software like HOMER
        """
        jobs=[]

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

            original_fastq1=fastq1 + "_original"
            original_fastq2=fastq2 + "_original"

            command = "mv {fastq1} {original_fastq1} && sed '/^@/s/\/[12]\>//g' {original_fastq1} > {fastq1}".format(fastq1 = fastq1, original_fastq1=original_fastq1)
            job_fastq1 = Job(input_files = [fastq1],
                    output_files = [fastq1],
                    name = "fastq1_readName_Edit." + readset.name,
                    command = command,
                    removable_files = [original_fastq1]
                    )

            command = "mv {fastq2} {original_fastq2} && sed '/^@/s/\/[12]\>//g' {original_fastq2} > {fastq2}".format(fastq2 = fastq2, original_fastq2=original_fastq2)
            job_fastq2 = Job(input_files = [fastq2],
                    output_files = [fastq2],
                    name = "fastq2_readName_Edit." + readset.name,
                    command = command,
                    removable_files = [original_fastq2]
                    )

            jobs.append(concat_jobs(job_fastq1, job_fastq2))

        return jobs


    def hicup_align(self):
        """
        Paired-end Hi-C reads are truncated, mapped and filtered using HiCUP. The resulting bam file is filtered for Hi-C artifacts and
        duplicated reads. It is ready for use as input for downstream analysis.

        For more detailed information about the HICUP process visit: [HiCUP] (https://www.bioinformatics.babraham.ac.uk/projects/hicup/overview/)
        """

        jobs = []

        # create hicup output directory:
        output_directory = "hicup_align"

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

            # create HiCUP configuration file:
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
            """.format(sample_output_dir = sample_output_dir,
                threads = config.param('hicup_align', 'threads'),
                Quiet = config.param('hicup_align', 'Quiet'),
                Keep = config.param('hicup_align', 'Keep'),
                Zip = config.param('hicup_align', 'Zip'),
                Bowtie2_path = os.path.expandvars(config.param('hicup_align', 'Bowtie2_path')),
                R_path = os.path.expandvars(config.param('hicup_align', 'R_path')),
                Genome_Index_hicup = os.path.expandvars(config.param('hicup_align', 'Genome_Index_hicup')),
                Genome_Digest = self.genome_digest,
                Format = config.param('hicup_align', 'Format'),
                Longest = config.param('hicup_align', 'Longest'),
                Shortest = config.param('hicup_align', 'Shortest'),
                fastq1 = fastq1,
                fastq2 = fastq2)

            ## write configFileContent to temporary file:
            fileName = "hicup_align." + readset.name + ".conf"
            with open(fileName, "w") as conf_file:
                conf_file.write(configFileContent)

            hicup_prefix = ".trim.pair1_2.hicup.bam" if config.param('hicup_align', 'Zip') == "1" else ".trim.pair1_2.hicup.sam"
            hicup_file_output = os.path.join("hicup_align", readset.name, readset.name + hicup_prefix)

            # hicup command
            ## delete directory if it exists since hicup will not run again unless old files are removed
            command = "rm -rf {sample_output_dir} && mkdir -p {sample_output_dir} && hicup -c {fileName}".format(sample_output_dir = sample_output_dir, fileName = fileName)
            job = Job(input_files = [fastq1, fastq2, fileName],
                    output_files = [hicup_file_output],
                    module_entries = [['hicup_align', 'module_bowtie2'], ['hicup_align', 'module_R'], ['hicup_align', 'module_mugqic_R_packages'], ['hicup_align', 'module_HiCUP']],
                    name = "hicup_align." + readset.name,
                    command = command,
                    removable_files = [fileName]
                    )

            jobs.append(job)

        return jobs




    def homer_tag_directory(self):
        """
        The bam file produced by HiCUP is used to create a tag directory using HOMER for further analysis that includes interaction matrix generation,
        compartments and identifying significant interactions.

        For more detailed information about the HOMER process visit: [HOMER] (http://homer.ucsd.edu/homer/interactions/index.html)
        """

        jobs = []

        output_directory = "homer_tag_directory"
        ## assuming that reads are trimmed and have a .trim. prefix, otherwise, edit code
        for readset in self.readsets:
            tagDirName = "_".join(("HTD", readset.name, self.enzyme))
            sample_output_dir = os.path.join(output_directory, tagDirName)
            hicup_prefix = ".trim.pair1_2.hicup.bam" if config.param('hicup_align', 'Zip') == "1" else ".trim.pair1_2.hicup.sam"
            hicup_file_output = os.path.join("hicup_align", readset.name, readset.name + hicup_prefix)
            
            ## homer command
            command_tagDir = "makeTagDirectory {sample_output_dir} {hicup_file_output},{hicup_file_output} -genome {genome} -restrictionSite {restriction_site} -checkGC".format(sample_output_dir = sample_output_dir, hicup_file_output = hicup_file_output, genome = config.param('DEFAULT', 'assembly'), restriction_site = self.restriction_site)

            command_archive = "mkdir -p {archive_output_dir} && mv -t {archive_output_dir} *random*.tsv *chrUn*.tsv chrM*.tsv chrY*.tsv".format(archive_output_dir = os.path.join(sample_output_dir, "archive"))

            #command_QcPlots = "makeTagDirectory {sample_output_dir} {hicup_file_output},{hicup_file_output} -genome {genome} -restrictionSite {restriction_site} -checkGC".format(sample_output_dir = sample_output_dir, hicup_file_output = hicup_file_output, genome = config.param('DEFAULT', 'assembly'), restriction_site = self.restriction_site)

            job = Job(input_files = [hicup_file_output],
                    output_files = [sample_output_dir],
                    module_entries = [["homer_tag_directory", "module_homer"]],
                    name = "homer_tag_directory." + readset.name,
                    command = command_tagDir + " && " + command_archive,
                    removable_files = []
                    )

            jobs.append(job)
        return jobs


    # def interaction_matrices(self):
    #     """
    #     IntraChromosomal interaction matrices, as well as genome wide interaction matrices are produced by Homer at resolutions defined in the ini config file
    #     For more detailed information about the HOMER matrices visit: [HOMER matrices] (http://homer.ucsd.edu/homer/interactions/HiCmatrices.html)
    #     """
        
    #     jobs = []

    #     chrs = config.param('interaction_matrices', 'chromosomes').split(",")
    #     res_chr = config.param('interaction_matrices', 'resolution_chr_matrix').split(",")
    #     res_genome = config.param('interaction_matrices', 'resolution_genome_matrix').split(",")
        
    #     output_directory = "interaction_matrices"
        
    #     for readset in self.readsets:
    #         tagDirName = "_".join(("HTD", readset.name, self.enzyme))
    #         homer_output_dir = os.path.join("homer_tag_directory", tagDirName)
    #         sample_output_dir = os.path.join(output_directory, readset.name)
            
    #         ## homer interaction matrices Chr command
    #         commandChr="makeTagDirectory {sample_output_dir} {hicup_file_output},{hicup_file_output} -genome {genome} -restrictionSite {restriction_site} -checkGC".format(sample_output_dir= sample_output_dir, hicup_file_output=hicup_file_output, genome=config.param('DEFAULT', 'assembly'), restriction_site= self.restriction_site)
    #         job = Job(input_files= [hicup_file_output],
    #                 output_files=[sample_output_dir],
    #                 module_entries= [["homer_tag_directory", "module_homer"]],
    #                 name= "homer_tag_directory." + readset.name,
    #                 command=command,
    #                 removable_files=[]
    #                 )

    #         jobs.append(job)
    #     return jobs


    def identify_compartments(self):
        """
        Genomic compartments are idetified using Homer at resolutions defined in the ini config file
        For more detailed information about the HOMER compartments visit: [HOMER compartments] (http://homer.ucsd.edu/homer/interactions/HiCpca.html)
        """
        pass

    def identify_TADs(self):
        """
        Topological associating Domains (TADs) are idetified using TopDom at resolutions defined in the ini config file
        For more detailed information about the TopDom visit: [TopDom] (https://www.ncbi.nlm.nih.gov/pubmed/26704975)
        """
        pass

    def identify_peaks(self):
        """
        Significant intraChromosomal interactions (peaks) are identified using Homer.
        For more detailed information about the Homer peaks visit: [Homer peaks] (http://homer.ucsd.edu/homer/interactions/HiCinteractions.html)
        """
        pass



    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.fastq_readName_Edit,
            self.hicup_align,
            self.homer_tag_directory
        ]

if __name__ == '__main__':
    HicSeq()
