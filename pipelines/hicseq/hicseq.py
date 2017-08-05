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


    @property
    def output_dirs(self):
        dirs = {'hicup_output_directory': 'hicup_align',
                'homer_output_directory': 'homer_tag_directory',
                'bams_output_directory': 'merged_bams',
                'matrices_output_directory': 'interaction_matrices',
                'cmpt_output_directory': 'identify_compartments',
                'TAD_output_directory': 'identify_TADs',
                'peaks_output_directory': 'identify_peaks'}
        return dirs



    def fastq_readName_Edit(self):
        """
        Removes the added /1 and /2 by picard's sam_to_fastq transformation to avoid issues with downstream software like HOMER
        """
        jobs=[]

        for readset in self.readsets:
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

            original_fastq1 = fastq1 + "original.gz"
            original_fastq2 = fastq2 + "original.gz"

            ## assumes reads in fastq file start with @; if not change
            command = "mv {fastq1} {original_fastq1} && zcat {original_fastq1} | sed '/^@/s/\/[12]\>//g' | gzip > {fastq1}".format(fastq1 = fastq1, original_fastq1=original_fastq1)
            job_fastq1 = Job(input_files = [fastq1],
                    output_files = [fastq1],
                    name = "fastq1_readName_Edit." + readset.name,
                    command = command,
                    removable_files = [original_fastq1]
                    )

            command = "mv {fastq2} {original_fastq2} && zcat {original_fastq2} | sed '/^@/s/\/[12]\>//g' | gzip > {fastq2}".format(fastq2 = fastq2, original_fastq2=original_fastq2)
            job_fastq2 = Job(input_files = [fastq2],
                    output_files = [fastq2],
                    name = "fastq2_readName_Edit." + readset.name,
                    command = command,
                    removable_files = [original_fastq2]
                    )

            jobs.extend([job_fastq1, job_fastq2])

        return jobs


    def hicup_align(self):
        """
        Paired-end Hi-C reads are truncated, mapped and filtered using HiCUP. The resulting bam file is filtered for Hi-C artifacts and
        duplicated reads. It is ready for use as input for downstream analysis.

        For more detailed information about the HICUP process visit: [HiCUP] (https://www.bioinformatics.babraham.ac.uk/projects/hicup/overview/)
        """

        jobs = []


        for readset in self.readsets:
            sample_output_dir = os.path.join(self.output_dirs['hicup_output_directory'], readset.name)
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
                    module_entries = [['hicup_align', 'module_bowtie2'], ['hicup_align', 'module_samtools'], ['hicup_align', 'module_R'], ['hicup_align', 'module_mugqic_R_packages'], ['hicup_align', 'module_HiCUP']],
                    name = "hicup_align." + readset.name,
                    command = command,
                    removable_files = [fileName]
                    )

            jobs.append(job)

        return jobs


    def samtools_merge_bams(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [samtools](http://samtools.sourceforge.net/).

        This step takes as input files the aligned bams/sams from the hicup_align step
        """

        jobs = []
        for sample in self.samples:
            sample_output = os.path.join(self.output_dirs['bams_output_directory'], sample.name + ".merged.bam")
            ## need to deal with sams if Zip=0
            readset_bams = [os.path.join(self.output_dirs['hicup_output_directory'], readset.name, readset.name + ".trim.pair1_2.hicup.bam") for readset in sample.readsets]

            mkdir_job = Job(command="mkdir -p " + self.output_dirs['bams_output_directory'])

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM.
            if len(sample.readsets) == 1:
                readset_bam = readset_bams[0]
                if os.path.isabs(readset_bam):
                    target_readset_bam = readset_bam
                else:
                    target_readset_bam = os.path.abspath(readset_bam)

                job = concat_jobs([
                    mkdir_job,
                    Job(input_files = readset_bams,
                    output_files = [sample_output],
                    command="ln -s -f " + target_readset_bam + " " + sample_output),
                ], name="symlink_readset_sample_bam." + sample.name)

            elif len(sample.readsets) > 1:

                samtools_merge = Job(input_files = readset_bams,
                    output_files = [sample_output],
                    module_entries = [['hicup_align', 'module_samtools']],
                    command = "samtools merge {sample_output} {input_bams}".format(sample_output = sample_output, input_bams = " ".join(map(str.strip, readset_bams)))
                    )

                job = concat_jobs([
                    mkdir_job,
                    samtools_merge
                ])
                job.name = "samtools_merge_bams." + sample.name

            jobs.append(job)

        return jobs



    def homer_tag_directory(self):
        """
        The bam file produced by HiCUP is used to create a tag directory using HOMER for further analysis that includes interaction matrix generation,
        compartments and identifying significant interactions.

        For more detailed information about the HOMER process visit: [HOMER] (http://homer.ucsd.edu/homer/interactions/index.html)
        """

        jobs = []

        ## assuming that reads are trimmed and have a .trim. prefix, otherwise, edit code
        for sample in self.samples:
            tagDirName = "_".join(("HTD", sample.name, self.enzyme))
            sample_output_dir = os.path.join(self.output_dirs['homer_output_directory'], tagDirName)
            hicup_file_output = os.path.join(self.output_dirs['bams_output_directory'], sample.name + ".merged.bam")
            archive_output_dir = os.path.join(sample_output_dir, "archive")
            QcPlots_output_dir = os.path.join(sample_output_dir, "HomerQcPlots")

            ## homer command
            command_tagDir = "rm -rf {sample_output_dir} && makeTagDirectory {sample_output_dir} {hicup_file_output},{hicup_file_output} -genome {genome} -restrictionSite {restriction_site} -checkGC".format(sample_output_dir = sample_output_dir, hicup_file_output = hicup_file_output, genome = config.param('DEFAULT', 'assembly'), restriction_site = self.restriction_site)

            command_QcPlots = "Rscript {script} {name} {workingDir} {outDir}".format(script = os.path.expandvars("${myhicseq}HomerHiCQcPlotGenerator.R"),  name = sample.name, workingDir = sample_output_dir, outDir = "HomerQcPlots")

            command_archive = "mkdir -p {archive_output_dir} && cd {sample_output_dir} && mv -t {archive} *random*.tsv *chrUn*.tsv *hap*.tsv chrM*.tsv chrY*.tsv".format(archive_output_dir = archive_output_dir, sample_output_dir = sample_output_dir, archive = "archive")


            job = Job(input_files = [hicup_file_output],
                    output_files = [sample_output_dir, archive_output_dir, QcPlots_output_dir],
                    module_entries = [["homer_tag_directory", "module_homer"], ["homer_tag_directory", "module_samtools"], ["homer_tag_directory", "module_R"]],
                    name = "homer_tag_directory." + sample.name,
                    command = command_tagDir + " && " + command_QcPlots + " && " + command_archive,
                    )

            jobs.append(job)
        return jobs



    def interaction_matrices_Chr(self):
        """
        IntraChromosomal interaction matrices are produced by Homer at resolutions defined in the ini config file and plotted by HiCPlotter.
        For more detailed information about the HOMER matrices visit: [HOMER matrices] (http://homer.ucsd.edu/homer/interactions/HiCmatrices.html)
        For more detailed information about HiCPlotter visit: [HiCPlotter] (https://github.com/kcakdemir/HiCPlotter)
        """

        jobs = []

        chrs = config.param('interaction_matrices_Chr', 'chromosomes').split(",")
        res_chr = config.param('interaction_matrices_Chr', 'resolution_chr').split(",")

        for sample in self.samples:
            tagDirName = "_".join(("HTD", sample.name, self.enzyme))
            homer_output_dir = os.path.join(self.output_dirs['homer_output_directory'], tagDirName)
            sample_output_dir_chr = os.path.join(self.output_dirs['matrices_output_directory'], sample.name, "chromosomeMatrices")

            # loop over chrs and res:
            for res in res_chr:
                for chr in chrs:

                    fileName = os.path.join(sample_output_dir_chr, "_".join((tagDirName, chr, res, "raw.txt")))
                    fileNameRN = os.path.join(sample_output_dir_chr, "_".join((tagDirName, chr, res, "rawRN.txt")))
                    commandChrMatrix="mkdir -p {sample_output_dir_chr} && analyzeHiC {homer_output_dir} -res {res} -raw -chr {chr} > {fileName} && cut -f 2- {fileName} > {fileNameRN}".format(sample_output_dir_chr = sample_output_dir_chr, homer_output_dir = homer_output_dir, res = res, chr = chr, fileName = fileName, fileNameRN = fileNameRN)
                    commandChrPlot = "HiCPlotter.py -f {fileNameRN} -n {name} -chr {chr} -r {res} -fh 0 -o {sample_output_dir_chr} -ptr 0 -hmc {hmc}".format(res=res, chr=chr, fileNameRN=fileNameRN, name=sample.name, sample_output_dir_chr=os.path.join(sample_output_dir_chr, "_".join((tagDirName, chr, res, "raw"))), hmc = config.param('interaction_matrices_Chr', 'hmc'))

                    job = Job(input_files = [homer_output_dir],
                            output_files = [fileNameRN, fileName],
                            module_entries = [["interaction_matrices_Chr", "module_homer"], ["interaction_matrices_Chr", "module_HiCPlotter"]],
                            name = "interaction_matrices_Chr." + sample.name + "_" + chr + "_res" + res,
                            command = commandChrMatrix + " && " + commandChrPlot,
                            removable_files = [fileName]
                            )

                    jobs.append(job)
        return jobs



    def interaction_matrices_genome(self):
        """
        Genomeside interaction matrices are produced by Homer at resolutions defined in the ini config file
        For more detailed information about the HOMER matrices visit: [HOMER matrices] (http://homer.ucsd.edu/homer/interactions/HiCmatrices.html)
        For more detailed information about HiCPlotter visit: [HiCPlotter] (https://github.com/kcakdemir/HiCPlotter)
        """

        jobs = []

        res_genome = config.param('interaction_matrices_genome', 'resolution_genome').split(",")


        for sample in self.samples:
            tagDirName = "_".join(("HTD", sample.name, self.enzyme))
            homer_output_dir = os.path.join(self.output_dirs['homer_output_directory'], tagDirName)
            sample_output_dir_genome = os.path.join(self.output_dirs['matrices_output_directory'], sample.name, "genomewideMatrices")

            for res in res_genome:
                fileName = os.path.join(sample_output_dir_genome, "_".join((tagDirName, "genomewide_Res", res,"raw.txt")))
                fileNameRN = os.path.join(sample_output_dir_genome, "_".join((tagDirName, "genomewide_Res", res,"rawRN.txt")))

                commandMatrix="mkdir -p {sample_output_dir_genome} && analyzeHiC {homer_output_dir} -res {res} -raw > {fileName} && cut -f 2- {fileName} > {fileNameRN}".format(sample_output_dir_genome = sample_output_dir_genome, homer_output_dir = homer_output_dir, res = res, chr = chr, fileName = fileName, fileNameRN = fileNameRN)

                commandPlot = "HiCPlotter.py -f {fileNameRN} -n {name} -chr Genome -r {res} -fh 0 -o {sample_output_dir_genome} -ptr 0 -hmc {hmc} -wg 1".format(res=res, chr=chr, fileNameRN=fileNameRN, name=sample.name, sample_output_dir_genome=os.path.join(sample_output_dir_genome, "_".join((tagDirName, "genomewide", "raw"))), hmc = config.param('interaction_matrices_Chr', 'hmc'))


                job = Job(input_files = [homer_output_dir],
                                output_files = [fileNameRN, fileName],
                                module_entries = [["interaction_matrices_genome", "module_homer"], ["interaction_matrices_genome", "module_HiCPlotter"]],
                                name = "interaction_matrices_genome." + sample.name  + "_res" + res,
                                command = commandMatrix + " && " + commandPlot,
                                removable_files = [fileName]
                                )

                jobs.append(job)
        return jobs



    def identify_compartments(self):
        """
        Genomic compartments are identified using Homer at resolutions defined in the ini config file
        For more detailed information about the HOMER compartments visit: [HOMER compartments] (http://homer.ucsd.edu/homer/interactions/HiCpca.html)
        """

        jobs = []

        res = config.param('identify_compartments', 'resolution_cmpt')

        for sample in self.samples:
            tagDirName = "_".join(("HTD", sample.name, self.enzyme))
            homer_output_dir = os.path.join(self.output_dirs['homer_output_directory'], tagDirName)
            sample_output_dir = os.path.join(self.output_dirs['cmpt_output_directory'], sample.name)
            fileName = os.path.join(sample_output_dir, sample.name + "_homerPCA_Res" + res)
            fileName_PC1 = os.path.join(sample_output_dir, sample.name + "_homerPCA_Res" + res + ".PC1.txt")
            fileName_Comp = os.path.join(sample_output_dir, sample.name + "_homerPCA_Res" + res + "_compartments")


            command = "mkdir -p {sample_output_dir} && runHiCpca.pl {fileName} {homer_output_dir} -res {res} -genome {genome}; findHiCCompartments.pl {fileName_PC1}  > {fileName_Comp}".format(sample_output_dir = sample_output_dir, fileName = fileName, homer_output_dir = homer_output_dir, res = res, genome = config.param('DEFAULT', 'assembly'), fileName_PC1 = fileName_PC1, fileName_Comp = fileName_Comp)

            job = Job(input_files = [homer_output_dir],
                    output_files = [fileName, fileName_PC1, fileName_Comp],
                    module_entries = [["identify_compartments", "module_homer"], ["identify_compartments", "module_R"]],
                    name = "identify_compartments." + sample.name + "_res" + res,
                    command = command,
                    removable_files = [fileName]
                    )

            jobs.append(job)
        return jobs



    def identify_TADs(self):
        """
        Topological associating Domains (TADs) are idetified using TopDom at resolutions defined in the ini config file
        For more detailed information about the TopDom visit: [TopDom] (https://www.ncbi.nlm.nih.gov/pubmed/26704975)
        """

        jobs = []

        chrs = config.param('identify_TADs', 'chromosomes').split(",")
        res_chr = config.param('identify_TADs', 'resolution_TADs').split(",")



        for sample in self.samples:
            sample_output_dir = os.path.join(self.output_dirs['TAD_output_directory'], sample.name)

            for res in res_chr:
                for chr in chrs:

                    input_matrix = os.path.join(self.output_dirs['matrices_output_directory'], sample.name, "chromosomeMatrices", "_".join(("HTD", sample.name, self.enzyme, chr, res, "rawRN.txt")))
                    tmp_matrix = input_matrix + ".MatA"
                    output_matrix = os.path.join(sample_output_dir, "_".join(("HTD", sample.name, self.enzyme, chr, res, "rawRN.MatA.TopDom")))

                    ## make TopDom R script:
                    FileContent = """
                    source("{script}")
                    TopDom(matrix.file="{tmp_matrix}", window.size={n}, outFile="{output_matrix}")
                    """.format(script = os.path.expandvars("${myhicseq}TopDom_v0.0.2.R"), tmp_matrix = tmp_matrix, n = config.param('identify_TADs', 'TopDom_n'), output_matrix = output_matrix)

                    ## write FileContent to temporary file:
                    fileName = "identify_TADs_TopDom." + sample.name + ".R"
                    with open(fileName, "w") as R_file:
                        R_file.write(FileContent)


                    command = "mkdir -p {sample_output_dir} && {script} {input} {res} && Rscript {Rscript}".format(sample_output_dir = sample_output_dir,script = os.path.expandvars("${myhicseq}CreateTopDomMat.sh"), input = input_matrix, res = res, Rscript = fileName)

                    job = Job(input_files = [input_matrix],
                            output_files = [output_matrix],
                            module_entries = [["identify_TADs", "module_R"]],
                            name = "identify_TADs." + sample.name + "_" + chr + "_res" + res,
                            command = command,
                            removable_files = [fileName, tmp_matrix]
                            )

                    jobs.append(job)
        return jobs



    def identify_peaks(self):
        """
        Significant intraChromosomal interactions (peaks) are identified using Homer.
        For more detailed information about the Homer peaks visit: [Homer peaks] (http://homer.ucsd.edu/homer/interactions/HiCinteractions.html)
        """

        jobs = []

        res = config.param('identify_peaks', 'resolution_pks')

        for readset in self.readsets:
            tagDirName = "_".join(("HTD", sample.name, self.enzyme))
            homer_output_dir = os.path.join(self.output_dirs['homer_output_directory'], tagDirName)
            sample_output_dir = os.path.join(self.output_dirs['peaks_output_directory'], sample.name)
            fileName = os.path.join(sample_output_dir, sample.name + "IntraChrInteractionsRes" + res + ".txt")
            fileName_anno = os.path.join(sample_output_dir, sample.name + "IntraChrInteractionsRes" + res + "_Annotated.txt")

            command = "mkdir -p {sample_output_dir} && findHiCInteractionsByChr.pl {homer_output_dir} -res {res} > {fileName} && annotateInteractions.pl {fileName} {genome} {fileName_anno}".format(sample_output_dir = sample_output_dir, homer_output_dir = homer_output_dir, res = res, fileName = fileName, genome = config.param('DEFAULT', 'assembly'), fileName_anno = fileName_anno)

            job = Job(input_files = [homer_output_dir],
                    output_files = [fileName, fileName_anno],
                    module_entries = [["identify_peaks", "module_homer"]],
                    name = "identify_peaks." + sample.name + "_res" + res,
                    command = command
                    )

            jobs.append(job)
        return jobs



    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.fastq_readName_Edit,
            self.hicup_align,
            self.samtools_merge_bams,
            self.homer_tag_directory,
            self.interaction_matrices_Chr,
            self.interaction_matrices_genome,
            self.identify_compartments,
            self.identify_TADs,
            self.identify_peaks
        ]

if __name__ == '__main__':
    HicSeq()
