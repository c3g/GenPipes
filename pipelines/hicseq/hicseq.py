#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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
import os
import re
import sys
import itertools

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from bfx import bash_cmd as bash
from core.config import config, SanitycheckError, _raise
from core.job import Job, concat_jobs, pipe_jobs
from pipelines import common
import utils.utils

from bfx import picard
from bfx import samtools
from bfx import hicup
from bfx import hicplotter
from bfx import homer
from bfx import multiqc
from bfx import genome
from bfx import bedtools
from bfx import chicago
from bfx import bedops
from bfx import tools
from bfx import topdom
from bfx import robustad
from bfx import hic
from bfx import hicrep
from bfx import quasar_qc

log = logging.getLogger(__name__)

class HicSeq(common.Illumina):
    """
    Hi-C Pipeline
    ==============

    Hi-C experiments allow researchers to understand chromosomal folding and structure using proximity ligation techniques.
    This pipeline analyzes both Hi-C experimental data (-t hic) and capture Hi-C data (-t capture).
    The Hi-C pipeline, selected using the "-t hic" parameter, starts by trimming adaptors and low quality bases.
    It then maps the reads to a reference genome using HiCUP.
    HiCUP first truncates chimeric reads, maps reads to the genome, then it filters Hi-C artifacts and removes read duplicates.
    Samples from different lanes are merged and a tag directory is created by Homer, which is also used to produce the interaction
    matrices and compartments. TopDom is used to predict topologically associating domains (TADs) and homer is used to identify
    significant interactions.

    The capture Hi-C pipeline, selected using the "-t capture" parameter, starts by trimming adaptors and low quality bases.
    It then maps the reads to a reference genome using HiCUP.
    HiCUP first truncates chimeric reads, maps reads to the genome, then it filters Hi-C artifacts and removes read duplicates.
    Samples from different lanes are merged and CHiCAGO is then used to filter capture-specific artifacts and call significant
    interactions. This pipeline also identifies enrichement of regulatory features when provided with ChIP-Seq marks. It can also
    return bed interctions with significant baits (bait_intersect step) or with captured interactions (capture_intersect step).

    An example of the Hi-C report for an analysis on public data (GM12878 Rao. et al.) is available for illustration purpose only:
    [Hi-C report](<url>).

    [Here](<url>) is more information about Hi-C pipeline that you may find interesting.
    """

    def __init__(self, protocol='hic'):
        self._protocol=protocol
        self.argparser.add_argument("-e", "--enzyme", help = "Restriction Enzyme used to generate Hi-C library (default DpnII)", choices = ["DpnII", "HindIII", "NcoI", "MboI", "Arima"], required=True, default="DpnII")
        self.argparser.add_argument("-t", "--type", help = "Hi-C experiment type (default hic)", choices = ["hic", "capture"], default="hic")
        super(HicSeq, self).__init__(protocol)

    @property
    def output_dirs(self):
        dirs = {'hicup_output_directory': 'alignment',
                'homer_output_directory': 'homer_tag_directory',
                'bams_output_directory': 'alignment',
                'matrices_output_directory': 'interaction_matrices',
                'cmpt_output_directory': 'compartments',
                'TAD_output_directory': 'TADs',
                'peaks_output_directory': 'peaks',
                'hicfiles_output_directory': 'hicFiles',
                'chicago_input_files': 'input_files',
                'chicago_output_directory': 'chicago',
                'intersect_ouput_directory': 'bed_intersect',
                'reproducible_score_output_directory': 'reproducibility_scores',
                'quality_score_output_directory': 'quality_scores'
                }
        return dirs

    @property
    def enzyme(self):
        return self.args.enzyme

    @property
    def restriction_site(self):
        """
        Sets the restriction enzyme recogntition site and genome digest location based on enzyme.
        """
        # Used only for Homer tag directory for QC of read location. For Arima, Homer doesn't accept multiple enzymes, use DpnII site
        if (self.enzyme == "DpnII") or (self.enzyme == "MboI") or (self.enzyme == "Arima"):
            restriction_site = "GATC"
        elif self.enzyme == "HindIII":
            restriction_site = "AAGCTT"
        elif self.enzyme == "NcoI":
            restriction_site = "CCATGG"
        else:
            _raise(SanitycheckError("Error: Selected Enzyme is not yet available for Hi-C analysis!"))
        return restriction_site

    @property
    def genome(self):
        genome_source = config.param('DEFAULT', 'source')
        if genome_source == "UCSC":
            genome = config.param('DEFAULT', 'assembly')
        else:
            genome = config.param('DEFAULT', 'assembly_synonyms')
        return genome

    @property
    def genome_digest(self):
        genome_digest = os.path.expandvars(config.param('hicup_align', "genome_digest_" + self.enzyme))
        return genome_digest

    def fastq_readName_Edit(self):
        """
        Removes the added /1 and /2 by picard's sam_to_fastq transformation to avoid issues with downstream software like HOMER.
        """
        jobs=[]

        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")

            if readset.run_type != "PAIRED_END":
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END for Hi-C analysis)!"))

            candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
            if readset.fastq1 and readset.fastq2:
                candidate_input_files.append([readset.fastq1, readset.fastq2])
            if readset.bam:
                candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam.strip())])
            [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            job_fastq1 = tools.sh_fastq_readname_edit(fastq1, self.output_dir, "fastq_readName_Edit.fq1." + readset.name)
            job_fastq1.samples = [readset.sample]

            job_fastq2 = tools.sh_fastq_readname_edit(fastq2, self.output_dir, "fastq_readName_Edit.fq2." + readset.name)
            job_fastq2.samples = [readset.sample]

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
            sample_output_dir = os.path.join(self.output_dirs['hicup_output_directory'], readset.sample.name, readset.name)
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")

            if readset.run_type != "PAIRED_END":
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END for Hi-C analysis)!"))

            candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
            if readset.fastq1 and readset.fastq2:
                candidate_input_files.append([readset.fastq1, readset.fastq2])
            if readset.bam:
                candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
            [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            # job_confFile = hicup.create_hicup_conf(readset.name, fastq1, fastq2, sample_output_dir, self.genome_digest)

            # job_hicup = hicup.hicup_run(readset.name, "hicup_align." + readset.name + ".conf", sample_output_dir, fastq1, fastq2, self.genome_digest)

            job = concat_jobs([
                    bash.mkdir(sample_output_dir),
                    hicup.hicup_run(
                                readset.name,
                                sample_output_dir,
                                fastq1,
                                fastq2,
                                self.genome_digest
                                )
                    ],
                    name = "hicup_align." + readset.name,
                    samples = [readset.sample]
                )

            # job = hicup.hicup_run(
            #                     readset.name,
            #                     sample_output_dir,
            #                     fastq1,
            #                     fastq2,
            #                     self.genome_digest
            #                     )

            jobs.append(job)

        return jobs

    def samtools_merge_bams(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [samtools](http://samtools.sourceforge.net/).

        This step takes as input files the aligned bams/sams from the hicup_align step.
        """

        jobs = []
        for sample in self.samples:
            sample_output = os.path.join(self.output_dirs['bams_output_directory'], sample.name, sample.name + ".merged.bam")
            readset_bams = [os.path.join(self.output_dirs['hicup_output_directory'], readset.sample.name, readset.name, readset.name + ".trim.pair1_2.fastq.gz.edited.hicup.bam") for readset in sample.readsets]

            mkdir_job = bash.mkdir(self.output_dirs['bams_output_directory'])

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM.
            if len(sample.readsets) == 1:
                readset_bam = readset_bams[0]
                if os.path.isabs(readset_bam):
                    target_readset_bam = readset_bam
                else:
                    target_readset_bam = os.path.relpath(readset_bam, os.path.dirname(sample_output))

                job = concat_jobs([
                    mkdir_job,
                    Job(
                        readset_bams,
                        [sample_output],
                        command="ln -s -f " + target_readset_bam + " " + sample_output
                    )
                ], name="symlink_readset_sample_bam." + sample.name)

            elif len(sample.readsets) > 1:

                samtools_merge = samtools.merge(sample_output, readset_bams)

                job = concat_jobs([
                    mkdir_job,
                    samtools_merge
                ])
                job.name = "samtools_merge_bams." + sample.name

            job.samples = [sample]
            jobs.append(job)

        return jobs

    def homer_tag_directory(self):
        """
        The bam file produced by HiCUP is used to create a tag directory using HOMER for further analysis that includes interaction matrix generation,
        compartments and identifying significant interactions.

        For more detailed information about the HOMER process visit: [HOMER] (http://homer.ucsd.edu/homer/interactions/index.html).
        """

        jobs = []

        makeDirTag_hic_other_options=config.param('homer_tag_directory', 'other_options', required=False)
        PEflag = eval(config.param('homer_tag_directory', 'illuminaPE'))

        for sample in self.samples:
            tagDirName = "_".join(("HTD", sample.name, self.enzyme))
            sample_output_dir = os.path.join(self.output_dirs['homer_output_directory'], tagDirName)
            hicup_file_output = os.path.join(self.output_dirs['bams_output_directory'], sample.name, sample.name + ".merged.bam")
            input_bam = ",".join([hicup_file_output, hicup_file_output])

            tagDir_job = homer.makeTagDir(sample_output_dir, input_bam, self.genome, self.restriction_site, PEflag, makeDirTag_hic_other_options)
            QcPlots_job = homer.hic_tagDirQcPlots(sample.name, sample_output_dir)
            archive_job = homer.archive_contigs(sample_output_dir)

            job = concat_jobs([tagDir_job, QcPlots_job, archive_job])
            job.name = "homer_tag_directory." + sample.name
            job.samples = [sample]

            jobs.append(job)
        return jobs

    def interaction_matrices_Chr(self):
        """
        IntraChromosomal interaction matrices are produced by Homer at resolutions defined in the ini config file and plotted by HiCPlotter.
        For more detailed information about the HOMER matrices visit: [HOMER matrices] (http://homer.ucsd.edu/homer/interactions/HiCmatrices.html).
        For more detailed information about HiCPlotter visit: [HiCPlotter] (https://github.com/kcakdemir/HiCPlotter).
        """

        jobs = []

        chrs = config.param('interaction_matrices_Chr', 'chromosomes')
        res_chr = config.param('interaction_matrices_Chr', 'resolution_chr').split(",")

        if chrs == "All":
            genome_dict = os.path.expandvars(config.param('DEFAULT', 'genome_dictionary', param_type='filepath'))
            chrs = genome.chr_names_conv(genome_dict)
        else:
            chrs = chrs.split(",")

        for sample in self.samples:
            tagDirName = "_".join(("HTD", sample.name, self.enzyme))
            homer_output_dir = os.path.join(self.output_dirs['homer_output_directory'], tagDirName)
            sample_output_dir_chr = os.path.join(self.output_dirs['matrices_output_directory'], sample.name, "chromosomeMatrices")

            # loop over chrs and res:
            for res in res_chr:
                for chr in chrs:

                    fileName = os.path.join(sample_output_dir_chr, "_".join((tagDirName, chr, res, "raw.txt")))
                    fileNameRN = os.path.join(sample_output_dir_chr, "_".join((tagDirName, chr, res, "rawRN.txt")))
                    fileNamePlot = os.path.join(sample_output_dir_chr, "".join((tagDirName,"_", chr,"_", res, "_raw-", chr, "\'.ofBins(0-\'*\')\'.", str(int(res)//1000), "K.jpeg")))
                    newFileNamePlot = os.path.join(sample_output_dir_chr, "".join((tagDirName,"_", chr,"_", res, "_raw-", chr, ".all.", str(int(res)//1000), "K.jpeg")))

                    jobMatrix = homer.hic_interactionMatrix_chr(sample.name, sample_output_dir_chr, homer_output_dir, res, chr, fileName, fileNameRN)
                    jobMatrix.samples = [sample]

                    jobPlot = hicplotter.intra_chrom_matrix_plot(
                        fileNameRN,
                        sample.name,
                        chr,
                        res,
                        os.path.join(sample_output_dir_chr, "_".join((tagDirName, chr, res, "raw"))),
                        fileNamePlot,
                        newFileNamePlot
                    )
                    jobPlot.samples = [sample]
                    jobs.extend([jobMatrix, jobPlot])

        return jobs

    def reproducibility_scores(self):

        """
        hic-rep is a R package for calculating the inter-chromosomal reproducibility score.
        Pairwise reproducibility scores for each chromosome pair in each sample pair are calculated using
        hic-rep at resolutions (bin size) defined in interaction_matrices_Chr step
        and other parameters defined in reproducibility_scores step of ini config file. All the scores are finally merged
        together and output a csv file with parameters used in analysis with chromosome number, reproducibility scores,
        standard deviation and smoothing value used for the analysis (in order to compare samples, smoothing value
        and the sequencing depth should be similar across samples).
        Down-sampling of samples can be performed using the down_sampling parameter in the ini config file.
        Correlation matrices and weight matrices can be saved using  corr=TRUE and weights=TRUE in ini config file
        for more information visit: [https://bioconductor.org/packages/release/bioc/html/hicrep.html]
        """

        jobs = []

        #Get defined chromosomes, resolution and other paramters from the config (.ini) file
        chrs = config.param('interaction_matrices_Chr', 'chromosomes')
        res_chr= [x.strip() for x in config.param('interaction_matrices_Chr', 'resolution_chr').split(",")]
        bound_width = config.param('reproducibility_scores', 'boundary_width')
        weights = config.param('reproducibility_scores', 'weights')
        down_sampling = config.param('reproducibility_scores', 'down_sampling')
        smooth = config.param('reproducibility_scores', 'h')
        corr = config.param('reproducibility_scores', 'corr')


        #Get the pairwise combinations of each sample (return an object list)
        pairwise_sample_combination = list(itertools.combinations(self.samples, 2))

        #run a loop for each resolution type and for each sample pair
        #compare each chromosome in the selected sample pair
        if chrs == "All":
            genome_dict = os.path.expandvars(config.param('DEFAULT', 'genome_dictionary', param_type='filepath'))
            chrs = genome.chr_names_conv(genome_dict)
        else:
            chrs = [x.strip() for x in chrs.split(",")]

        input_files_for_merging=[]
        tsv_files_for_merging=[]

        #If there is only one sample, an step will skip by giving a warning
        if len(self.samples) > 1:

            output_dir = self.output_dirs['reproducible_score_output_directory']
            temp_dir = "temp"
            for res in res_chr:
                for sample in pairwise_sample_combination:
                    #create a job list inorder to create a concat jobs (create one job for all chrms for one sample pair
                    #comparison
                    job_all_chr = []
                    for chromosome in chrs:

                        #get sample names with relative file path
                        #since pairwise_sample_combination is a list with two elements, in each index, 0 and 1 have been
                        #used to call the each sample
                        input_sample1_file_path = os.path.join(
                            self.output_dirs['matrices_output_directory'],
                            sample[0].name,
                            "chromosomeMatrices",
                            "_".join(("HTD", sample[0].name, self.enzyme, chromosome, res, "rawRN.txt"))
                            )

                        input_sample2_file_path = os.path.join(self.output_dirs['matrices_output_directory'],
                            sample[1].name,
                            "chromosomeMatrices",
                            "_".join(("HTD", sample[1].name, self.enzyme, chromosome, res, "rawRN.txt"))
                            )
                        out_dir = os.path.join(
                            output_dir,
                            temp_dir,
                           "_".join((sample[0].name, "vs",  sample[1].name))
                           )

                        output_file = os.path.join(out_dir, "hicrep_" + sample[0].name + "_vs_" + sample[1].name + "_" + chromosome + "_" + res + "_res_" + smooth + "_" + bound_width + "_" + down_sampling + ".tmp")

                        job_chr = concat_jobs([
                            bash.mkdir(out_dir),
                            hicrep.calculate_reproducible_score(
                                output_file,
                                sample[0].name,
                                sample[1].name,
                                input_sample1_file_path,
                                input_sample2_file_path,
                                chromosome,
                                res,
                                bound_width,
                                weights,
                                corr,
                                down_sampling,
                                smooth
                                )
                            ])

                        job_chr.samples = sample

                        input_files_for_merging.append(output_file)
                        job_all_chr.append(job_chr)
                        job = concat_jobs(job_all_chr)
                        job.name = "_".join(("reproducibility_scores.hicrep", sample[0].name, "vs", sample[1].name, res))
                        job.removable_files = [output_file]

                    jobs.append(job)
                    tsv_output = os.path.join(
                        output_dir,
                        temp_dir,
                        ".".join(("_".join((sample[0].name,"vs",sample[1].name, "res", res, smooth, bound_width, down_sampling)), "tsv"))
                        )
                    #Storing output files for next step in a list
                    tsv_files_for_merging.append(tsv_output)
                #create a job for merginh .tsv file with individual reproducibnility score for each pairwise comparison
                job_merge = hicrep.merge_tmp_files(
                                                    input_files_for_merging,
                                                    tsv_files_for_merging,
                                                    output_dir,
                                                    res,
                                                    smooth,
                                                    bound_width,
                                                    down_sampling,
                                                    temp_dir
                                                    )
                job_merge.samples = self.samples
                job_merge.name = "merge_hicrep_scores." + res
                jobs.append(job_merge)

            #finally all the .tsv files are merged and create a one file (.csv) with all the values
            output_csv_file = "hicrep_combined_reproducibility_scores.csv"
            job_tsv_merge = hicrep.merge_tsv(
                tsv_files_for_merging,
                output_dir,
                output_csv_file,
                temp_dir
                )
            job_tsv_merge.name = "merge_hicrep_scores"
            job_tsv_merge.samples = self.samples
            job_tsv_merge.removable_files = [os.path.join(output_dir, temp_dir)]
            jobs.append(job_tsv_merge)
        else:
            log.info("You need to provide at least two samples to generate a reproducible score... skipping")

        return jobs

    def quality_scores(self):

        """
        Quality score per chromosome for each sample is calculated using QUASAR-QC at all resolutions
        and sequencing depths (coverages) and down_sampling value (coverage) defined in quality_scores step of ini config file
        QUASAR-QC is a part of the hifive hic-seq analysis suite
        for more information visit: [http://hifive.docs.taylorlab.org/en/latest/quasar_scoring.html]
        """
        jobs = []

        # Get defined chromosomes and resolution from the config (.ini) file
        chrs = config.param('interaction_matrices_Chr', 'chromosomes')
        res_chr= [x.strip() for x in config.param('interaction_matrices_Chr', 'resolution_chr').split(",")]
        genome_fasta = config.param('DEFAULT', 'genome_fasta')
        chrom_lengths = ".".join((genome_fasta, "fai"))
        output_dir = self.output_dirs['quality_score_output_directory']
        temp_dir="temp"
        output_fend_file=""
        for res in res_chr:
            #First of all, a fend object should be created and this can be used with all the other samples.
            # Create fend object using chromosome lengths. only input file is chromosome lengths for respective genome
            output_fend_file = os.path.join(output_dir, temp_dir, ".".join(("binned_by_chr_length", res, "fends")))
            job_fend = quasar_qc.create_fend_object(chrom_lengths, output_fend_file, output_dir, res, temp_dir)
            job_fend.samples = self.samples
            job_fend.name = "quality_scores_quasar_fend"
            jobs.append(job_fend)

        if chrs == "All":
            genome_dict = os.path.expandvars(config.param('DEFAULT', 'genome_dictionary', param_type='filepath'))
            chrs = genome.chr_names_conv(genome_dict)
        else:
            chrs = [x.strip() for x in chrs.split(",")]

        #Restrutcturing interaction matrices suitable for QUASAR-QC input
        #First need to add a lable as (chr_number:bin_start - bin_end) and set it as the row name (first column)
        #Remove the exisiting column headers and first 3 columns
        #loop through all the samples and chromsomes
        #these files will be used in next steps
        quasr_temp_files = []
        quasr_temp_files.append(output_fend_file)
        for sample in self.samples:
            job_restructure = []
            for res in res_chr:
                for chromosome in chrs:

                    input_file_path = os.path.join(self.output_dirs['matrices_output_directory'],
                                                          sample.name,
                                                          "chromosomeMatrices",
                                                          "_".join((
                                                              "HTD", sample.name, self.enzyme, chromosome, res,
                                                              "rawRN.txt")))
                    output_file_path = os.path.join(output_dir, temp_dir, "_".join((
                                                              "HTD", sample.name, self.enzyme, chromosome, res,
                                                              "rawRN.txt.temp")))
                    quasr_temp_files.append(output_file_path)
                    job_matrix_restucture = quasar_qc.restructure_matrix( input_file_path, output_file_path,  output_dir, res, temp_dir)

                    job_matrix_restucture.samples = [sample]
                    job_restructure.append(job_matrix_restucture)
                    job = concat_jobs(job_restructure)
                    job.name = "_".join(("quality_scores.quasar", "restructuring",sample.name, res, self.enzyme))

            jobs.append(job)


        #Perform the QUASAR-QC analysis using the generated temp files
        #loop through each sample with same enzyme and one or many resolutions user has specified.
        #user can specify the downsampling value too. If needed the full coverage (all reads), use 0 as the downsampling value

        quasr_res = config.param('quality_scores', 'resolution_chr')
        quasr_coverage = config.param('quality_scores', 'down_sampling')

        quasr_report_files = []
        quasr_prefix = "quasarqc_res"

        for sample in self.samples:
            for res in res_chr:
                input_files = os.path.join(output_dir, temp_dir, "_".join((
                                                          "HTD", sample.name, self.enzyme, "*", res,
                                                          "rawRN.txt.temp")))

                output_file_path = os.path.join(output_dir, temp_dir, "_".join((
                    sample.name, quasr_prefix ,res, self.enzyme)))
                output_report_file = "_".join((output_file_path, "report.txt"))
                quasr_report_files.append(output_report_file)
                job_qusarqc = quasar_qc.quality_analysis(quasr_temp_files, input_files, output_file_path, output_dir, output_fend_file, quasr_res, quasr_coverage, output_report_file)

                job_qusarqc.samples = [sample]
                job_qusarqc.name = "_".join(("quality_scores.quasarqc",sample.name, res, self.enzyme))
                jobs.append(job_qusarqc)

        for res in res_chr:
            for q_res in quasr_res.split(","):
                q_res = int(q_res)/1000
                q_res = "".join((str(q_res) , "Kb" ))
                input_files = [s for s in quasr_report_files if "_".join((quasr_prefix, res, self.enzyme)) in s]
                output_file = ".".join((os.path.join(output_dir, "_".join((
                    "quasar_qc_final_report_res", res, "quasarqc_res", q_res, self.enzyme))), "tsv"))
                job_qusarqc_final_report = quasar_qc.merge_all_reports(output_dir, input_files, output_file, res, self.enzyme, quasr_prefix, temp_dir, q_res)
                job_qusarqc_final_report.samples = self.samples
                job_qusarqc_final_report.name ="_".join(("quality_scores_merge.quasarqc_final_report","res", res, "quasar_res", q_res))
                job_qusarqc_final_report.removable_files = [os.path.join(output_dir,temp_dir)]
                jobs.append(job_qusarqc_final_report)

        return jobs


    def interaction_matrices_genome(self):
        """
        Genomewide interaction matrices are produced by Homer at resolutions defined in the ini config file.
        For more detailed information about the HOMER matrices visit: [HOMER matrices] (http://homer.ucsd.edu/homer/interactions/HiCmatrices.html).
        For more detailed information about HiCPlotter visit: [HiCPlotter] (https://github.com/kcakdemir/HiCPlotter).
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

                jobMatrix = homer.hic_interactionMatrix_genome (sample.name, sample_output_dir_genome, homer_output_dir, res, fileName, fileNameRN)
                jobMatrix.samples = [sample]

                jobPlot = hicplotter.genome_wide_matrix_plot(
                    fileNameRN,
                    sample.name,
                    res,
                    os.path.join(sample_output_dir_genome, "_".join((tagDirName, "genomewide", "raw"))),
                    os.path.join(sample_output_dir_genome, tagDirName + "_genomewide_raw-WholeGenome-" + str(int(res)//1000) + "K.jpeg")
                )
                jobPlot.samples = [sample]

                jobs.extend([jobMatrix, jobPlot])

        return jobs

    def identify_compartments(self):
        """
        Genomic compartments are identified using Homer at resolutions defined in the ini config file.
        For more detailed information about the HOMER compartments visit: [HOMER compartments] (http://homer.ucsd.edu/homer/interactions/HiCpca.html).
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

            job = homer.hic_compartments (sample.name, sample_output_dir, fileName, homer_output_dir, res, self.genome, fileName_PC1, fileName_Comp, 3)
            job.samples = [sample]

            jobs.append(job)

        return jobs

    def identify_TADs_TopDom(self):
        """
        Topological associating Domains (TADs) are identified using TopDom at resolutions defined in the ini config file.
        For more detailed information about the TopDom visit: [TopDom] (https://www.ncbi.nlm.nih.gov/pubmed/26704975).
        """

        jobs = []

        chrs = config.param('identify_TADs', 'chromosomes')
        res_chr = config.param('identify_TADs', 'resolution_TADs').split(",")

        if chrs == "All":
            genome_dict = os.path.expandvars(config.param('DEFAULT', 'genome_dictionary', param_type='filepath'))
            chrs = genome.chr_names_conv(genome_dict)
        else:
            chrs = chrs.split(",")

        for sample in self.samples:
            sample_output_dir = os.path.join(self.output_dirs['TAD_output_directory'], sample.name, "TopDom")

            for res in res_chr:
                for chr in chrs:

                    input_matrix = os.path.join(self.output_dirs['matrices_output_directory'], sample.name, "chromosomeMatrices", "_".join(("HTD", sample.name, self.enzyme, chr, res, "rawRN.txt")))
                    tmp_matrix = input_matrix + ".MatA"
                    output_matrix = os.path.join(sample_output_dir, "_".join(("HTD", sample.name, self.enzyme, chr, res, "rawRN.MatA.TopDom")))
                    output_script = "identify_TADs_TopDom." + sample.name + "_" + chr + "_res" + res + ".R"

                    job_inputFile = concat_jobs(
                        [
                            bash.mkdir(sample_output_dir),
                            topdom.create_input(input_matrix, tmp_matrix, output_matrix, output_script, res)
                        ],
                        name="identify_TADs.TopDom_create_input." + sample.name + "_" + chr + "_res" + res,
                        samples=[sample]
                    )

                    job_TADs = topdom.call_TADs(tmp_matrix, output_matrix, output_script)
                    job_TADs.name = "identify_TADs.TopDom_call_TADs." + sample.name + "_" + chr + "_res" + res
                    job_TADs.samples = [sample]

                    jobs.extend([
                        job_inputFile,
                        job_TADs
                    ])

        return jobs

    def identify_TADs_RobusTAD(self):
        """
        Topological associating Domain (TAD) scores are calculated using RobusTAD for every bin in the genome.
        RobusTAD is resolution-independant and will use the first resolution in "resolution_TADs"  under [identify_TADs] in the ini file.
        For more detailed information about the RobusTAD visit: [RobusTAD] (https://github.com/rdali/RobusTAD).
        """

        jobs = []

        chrs = config.param('identify_TADs', 'chromosomes')
        res = config.param('identify_TADs', 'resolution_TADs').split(",")[0]

        if chrs == "All":
            genome_dict = os.path.expandvars(config.param('DEFAULT', 'genome_dictionary', param_type='filepath'))
            chrs = genome.chr_names_conv(genome_dict)
        else:
            chrs = chrs.split(",")

        for sample in self.samples:
            sample_output_dir = os.path.join(self.output_dirs['TAD_output_directory'], sample.name, "RobusTAD")

            for chr in chrs:

                input_matrix = os.path.join(self.output_dirs['matrices_output_directory'], sample.name, "chromosomeMatrices", "_".join(("HTD", sample.name, self.enzyme, chr, res, "rawRN.txt")))
                prefix = os.path.splitext(os.path.basename(input_matrix))[0]
                output_Scores = os.path.join(sample_output_dir, "".join(("BoundaryScores_", prefix, "_binSize" , str(int(res)//1000) ,"_minW250_maxW500_minRatio1.5.txt")))
                output_calls = os.path.join(sample_output_dir, "".join(("TADBoundaryCalls_", prefix, "_binSize" , str(int(res)//1000) ,"_minW250_maxW500_minRatio1.5_threshold0.2.txt")))

                job = concat_jobs([
                    bash.mkdir(sample_output_dir),
                    robustad.call_TADs(input_matrix, sample_output_dir, res)
                ])
                job.name = "identify_TADs.RobusTAD." + sample.name + "_" + chr
                job.samples = [sample]

                jobs.append(job)

        return jobs

    def identify_peaks(self):
        """
        Significant intraChromosomal interactions (peaks) are identified using Homer.
        For more detailed information about the Homer peaks visit: [Homer peaks] (http://homer.ucsd.edu/homer/interactions/HiCinteractions.html).
        """

        jobs = []

        res = config.param('identify_peaks', 'resolution_pks')

        for sample in self.samples:
            tagDirName = "_".join(("HTD", sample.name, self.enzyme))
            homer_output_dir = os.path.join(self.output_dirs['homer_output_directory'], tagDirName)
            sample_output_dir = os.path.join(self.output_dirs['peaks_output_directory'], sample.name)
            fileName = os.path.join(sample_output_dir, sample.name + "IntraChrInteractionsRes" + res + ".txt")
            fileName_anno = os.path.join(sample_output_dir, sample.name + "IntraChrInteractionsRes" + res + "_Annotated")

            job = homer.hic_peaks(sample.name, sample_output_dir, homer_output_dir, res, self.genome, fileName, fileName_anno, 3)
            job.samples = [sample]

            jobs.append(job)

        return jobs

    def create_hic_file(self):
        """
        A .hic file is created per sample in order to visualize in JuiceBox, WashU epigenome browser or as input for other tools.
        For more detailed information about the JuiceBox visit: [JuiceBox] (http://www.aidenlab.org/software.html).
        """

        jobs = []

        for sample in self.samples:
            sample_input = os.path.join(self.output_dirs['bams_output_directory'], sample.name, sample.name + ".merged.bam")
            sortedBamPrefix = re.sub("\.merged.bam", ".merged.sorted", sample_input.strip())
            sortedBam = sortedBamPrefix + ".bam"
            hic_output = os.path.join(self.output_dirs['hicfiles_output_directory'], sample.name + ".hic")

            job = concat_jobs([
                bash.mkdir(self.output_dirs['hicfiles_output_directory']),
                samtools.sort(sample_input, sortedBamPrefix, sort_by_name=True),
                hic.create_input(sortedBam, sample.name),
                hic.create_hic(sample.name + ".juicebox.input.sorted", hic_output, self.genome)
            ])
            job.name = "create_hic_file." + sample.name
            job.samples = [sample]

            jobs.append(job)

        return jobs

    def multiqc_report(self):
        """
        A quality control report for all samples is generated.
        For more detailed information about the MultiQc visit: [MultiQc] (http://multiqc.info/).
        """
        ## set multiQc config file so we can customize one for every pipeline:

        jobs = []

        yamlFile = os.path.expandvars(config.param('multiqc_report', 'MULTIQC_CONFIG_PATH'))
        input_files = [os.path.join(self.output_dirs['bams_output_directory'], sample.name, sample.name + ".merged.bam") for sample in self.samples]
        job = multiqc.mutliqc_run(yamlFile, input_files)
        job.samples = self.samples

        jobs.append(job)

        return jobs

        ## capture HiC methods:

    def create_rmap_file(self):
        """
        rmap file for Chicago capture analysis is created using the hicup digestion file.
        """
        ## return 1 rmap per enzyme

        output = os.path.join(self.output_dirs['chicago_input_files'], self.enzyme + ".Initialrmap")
        sorted_output = re.sub("\.Initialrmap", ".sorted.rmap", output)
        input_file = self.genome_digest

        job = concat_jobs([
            bash.mkdir(self.output_dirs['chicago_input_files']),
            tools.sh_create_rmap(input_file, output, "create_rmap_file." + self.enzyme),
            bedops.sort_bed(output, sorted_output)
        ])
        job.name = "create_rmap_file." + os.path.basename(input_file)
        job.samples = self.samples

        return [job]

    def create_baitmap_file(self):
        """
        baitmap file for Chicago capture analysis is created using the created rmap file and the probe capture bed file.
        """
        ## return 1 baitmap per enzyme/capture array combination

        input_rmap = os.path.join(self.output_dirs['chicago_input_files'], self.enzyme + ".sorted.rmap")
        input_bait = config.param('create_baitmap_file', "baitBed")
        sorted_input_bait = re.sub("\.bed", ".sorted.bed", input_bait)
        output_file_name = re.sub("\.bed", "", os.path.basename(input_bait)) + "_" + self.enzyme + ".baitmap"
        output_file = os.path.join(self.output_dirs['chicago_input_files'], output_file_name)
        annotation = config.param('create_baitmap_file', "annotation")

        job = concat_jobs([
            bedops.sort_bed(input_bait, sorted_input_bait),
            bedtools.intersect_beds(input_rmap, sorted_input_bait, output_file + ".tmp", "-wa -u"),
            tools.sh_create_baitmap(input_bait, sorted_input_bait, annotation, output_file)
        ])
        job.name = "create_baitmap_file." + output_file_name
        job.samples = self.samples

        return [job]

    def create_design_files(self):
        """
        design files (nperbin file (.npb), nbaitsperbin file (.nbpb), proxOE file (.poe)) for Chicago capture analysis are created using the rmap file and the baitmap file.
        """

        rmapfile = os.path.join(self.output_dirs['chicago_input_files'], self.enzyme + ".sorted.rmap")
        baitmapfile = os.path.join(self.output_dirs['chicago_input_files'], os.path.basename(re.sub("\.bed$", "", config.param('create_baitmap_file', "baitBed")) + "_" + self.enzyme + ".baitmap"))
        other_options = config.param('create_design_files', 'other_options', required = False)
        designDir = self.output_dirs['chicago_input_files']
        outfilePrefix = os.path.join(self.output_dirs['chicago_input_files'], os.path.basename(re.sub("\.bed$", "", config.param('create_baitmap_file', "baitBed")) + "_" + self.enzyme))
        job = chicago.makeDesignFiles(rmapfile, baitmapfile, outfilePrefix, designDir, other_options)
        job.samples = self.samples

        return [job]

    def create_input_files(self):
        """
        input file (sample.chinput) for Chicago capture analysis is created using the rmap file, the baitmap file and the hicup aligned bam.
        """

        jobs = []
        rmapfile = os.path.join(self.output_dirs['chicago_input_files'], self.enzyme + ".sorted.rmap")
        baitmapfile = os.path.join(self.output_dirs['chicago_input_files'], os.path.basename(re.sub("\.bed$", "", config.param('create_baitmap_file', "baitBed")) + "_" + self.enzyme + ".baitmap"))
        other_options = config.param('create_input_files', 'other_options', required = False)

        for sample in self.samples:
            name = os.path.join(self.output_dirs['chicago_input_files'], sample.name)
            bam = os.path.join(self.output_dirs['bams_output_directory'], sample.name, sample.name + ".merged.bam")
            job = chicago.bam2chicago(bam, baitmapfile, rmapfile, name, other_options)
            job.samples = [sample]

            jobs.append(job)

        return jobs

    def runChicago(self):
        """
        Chicago is run on capture data. Chicago will filter capture hic artifacts and identify significant interactions. It will output data as a bed file and will also output SeqMonk and WashU tracks.
        For more detailed information about the Chicago, including how to interpret the plots, please visit: [Chicago](https://bioconductor.org/packages/release/bioc/vignettes/Chicago/inst/doc/Chicago.html).
        """

        jobs = []
        design_dir = self.output_dirs['chicago_input_files']
        output_dir = self.output_dirs['chicago_output_directory']
        design_file_prefix = os.path.basename(re.sub("\.bed$", "", config.param('create_baitmap_file', "baitBed"))
                                              + "_" + self.enzyme)
        other_options = config.param('runChicago', 'other_options', required=False)

        for sample in self.samples:
            job = chicago.runChicago(design_dir, sample.name, output_dir, design_file_prefix, other_options)
            job.samples = [sample]

            jobs.append(job)

        return jobs

    def runChicago_featureOverlap(self):
        """
        Runs the feature enrichement of chicago significant interactions.
        For more detailed information about the Chicago, including how to interpret the plots, please visit: [Chicago](https://bioconductor.org/packages/release/bioc/vignettes/Chicago/inst/doc/Chicago.html).
        """

        jobs = []
        design_dir = self.output_dirs['chicago_input_files']
        output_dir = self.output_dirs['chicago_output_directory']
        design_file_prefix = os.path.basename(re.sub("\.bed$", "", config.param('create_baitmap_file', "baitBed")) + "_" + self.enzyme)
        other_options = config.param('runChicago_featureOverlap', 'other_options', required = False)
        features_file = config.param('runChicago_featureOverlap', 'features_file', required = True)

        if features_file != "None":
            for sample in self.samples:
                job = chicago.runChicago_featureOverlap(features_file, sample.name, output_dir, design_file_prefix, other_options)
                job.samples = [sample]

                jobs.append(job)

        return jobs

    def bait_intersect(self):
        """
        provided with a bed file, for example a bed of GWAS snps or features of interest, this method returns the lines in the bed file that intersect with the baits that have significant interactions.
        Input bed must have 4 columns (<chr> <start> <end> <annotation>) and must be tab separated.
        """

        jobs = []
        chicago_output_dir = self.output_dirs['chicago_output_directory']
        intersect_output_dir = self.output_dirs['intersect_ouput_directory']
        other_options = config.param('bait_intersect', 'other_options', required = False)
        features_file = config.param('bait_intersect', 'features_file', required = True)

        sorted_features_file = os.path.splitext(features_file)[0] + ".sorted.bed"
        output_dir = os.path.join(chicago_output_dir, intersect_output_dir)

        if features_file != "None":

            job_sort_features_file = bedops.sort_bed(features_file, sorted_features_file)

            job = concat_jobs([
                bash.mkdir(output_dir),
                job_sort_features_file
            ])
            job.name = "bait_intersect.sort_feature." + os.path.basename(sorted_features_file)
            job.removable_files = [sorted_features_file]
            job.samples = self.samples

            jobs.append(job)

            for sample in self.samples:
                sample_output_dir = os.path.join(chicago_output_dir, sample.name, "data")
                ibed_file = os.path.join(sample_output_dir, sample.name + ".ibed")
                sorted_ibed_file = re.sub("\.ibed$", ".bait.sorted.bed", ibed_file)

                output_file_prefix = os.path.join(output_dir, os.path.splitext(os.path.basename(features_file))[0] + "_" + os.path.splitext(os.path.basename(ibed_file))[0])

                job_extract_bait_bed = tools.sh_extract_bait_bed(ibed_file, sample.name)
                job_sort_ibed = bedops.sort_bed(ibed_file + ".bait", sorted_ibed_file)
                job_intersect = bedtools.intersect_beds(sorted_features_file, sorted_ibed_file, output_file_prefix + ".tmp", "-wa -u")
                job_bedopsMap = bedops.bedmap_echoMapId(output_file_prefix + ".tmp", sorted_ibed_file, output_file_prefix + ".bait_intersect.bed")

                job = concat_jobs([
                    job_extract_bait_bed,
                    job_sort_ibed,
                    job_intersect,
                    job_bedopsMap
                ])
                job.name = "bait_intersect." + sample.name
                job.removable_files = [sorted_ibed_file, sorted_features_file, output_file_prefix + ".tmp", ibed_file + ".bait"]
                job.samples = [sample]

                jobs.append(job)

        return jobs

    def capture_intersect(self):
        """
        provided with a bed file, for example a bed of GWAS snps or features of interest, this method returns the lines in the bed file that intersect with the captured ends ("Other Ends") that have significant interactions.
        Input bed must have 4 columns (<chr> <start> <end> <annotation>) and must be tab separated.
        """

        jobs = []
        chicago_output_dir = self.output_dirs['chicago_output_directory']
        intersect_output_dir = self.output_dirs['intersect_ouput_directory']
        other_options = config.param('capture_intersect', 'other_options', required = False)
        features_file = config.param('capture_intersect', 'features_file', required = True)

        sorted_features_file = os.path.splitext(features_file)[0] + ".sorted.bed"
        output_dir = os.path.join(chicago_output_dir, intersect_output_dir)

        if features_file != "None":
            job_sort_features_file = bedops.sort_bed(features_file, sorted_features_file)

            job = concat_jobs([bash.mkdir(output_dir), job_sort_features_file])
            job.name = "capture_intersect.sort_feature." + os.path.basename(sorted_features_file)
            job.removable_files = [sorted_features_file]
            job.samples = self.samples

            jobs.append(job)

            for sample in self.samples:
                sample_output_dir = os.path.join(chicago_output_dir, sample.name, "data")
                ibed_file = os.path.join(sample_output_dir, sample.name + ".ibed")
                sorted_ibed_file = re.sub("\.ibed$", ".capture.sorted.bed", ibed_file)

                output_file_prefix = os.path.join(output_dir, os.path.splitext(os.path.basename(features_file))[0] + "_" + os.path.splitext(os.path.basename(ibed_file))[0])

                job_extract_capture_bed = tools.sh_extract_capture_bed(ibed_file, sample.name)
                job_sort_ibed = bedops.sort_bed(ibed_file + ".capture", sorted_ibed_file)
                job_intersect = bedtools.intersect_beds(sorted_features_file, sorted_ibed_file, output_file_prefix + ".tmp", "-wa -u")
                job_bedopsMap = bedops.bedmap_echoMapId(output_file_prefix + ".tmp", sorted_ibed_file, output_file_prefix + ".capture_intersect.bed")

                job = concat_jobs([
                    job_extract_capture_bed,
                    job_sort_ibed,
                    job_intersect,
                    job_bedopsMap
                ])
                job.name = "capture_intersect." + sample.name
                job.removable_files = [sorted_ibed_file, output_file_prefix + ".tmp", ibed_file + ".capture"]
                job.samples = [sample]

                jobs.append(job)

        return jobs

    @property
    def steps(self):
        return [
            [self.samtools_bam_sort,
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
            self.identify_TADs_TopDom,
            self.identify_TADs_RobusTAD,
            self.identify_peaks,
            self.create_hic_file,
            self.reproducibility_scores,
            self.quality_scores,
            self.cram_output,
            self.multiqc_report
            ],
            [self.samtools_bam_sort,
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
            self.runChicago,
            self.runChicago_featureOverlap,
            self.bait_intersect,
            self.capture_intersect,
            self.create_hic_file,
            self.multiqc_report,
            self.cram_output
            ]
        ]

if __name__ == '__main__':
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        HicSeq(protocol=['hic','capture'])
