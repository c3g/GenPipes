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
import argparse
import logging
import os
import sys
import csv
from operator import attrgetter


# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import global_config_parser, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs

from pipelines import common

from bfx import bigwiginfo
from bfx import chromimpute
from bfx import epigeec
from bfx import epiqc_reports
from bfx.readset import parse_illumina_readset_file
import utils.utils
from shutil import copyfile
from bfx import genome

#from pipelines.chipseq import chipseq
from bfx import bash_cmd as bash
log = logging.getLogger(__name__)


class EpiQC(common.Illumina):
    """
        EpiQC Pipeline
        ==============

        EpiQC is a quality control pipeline for signal files (bigwig) generated from ChIP-Seq. The pipeline does a series of calculations on
        these files to assess the quality of ChIp-Seq data. Four different metrics are computed from a single bigwig file.
        First, BigWigInfo will be used to initial quality check on signal tracks
        Second, ChromImpute imputes signal tracks for the given chromosome (currently only chr1 is supported) and using these imputed files, EpiQC computes 2 other metrics.
        Third, SignalToNoise measurement by calculating the proportion of signal in top bins
        And Fourth, the pipeline creates a heatmap from the correlation matrix obtained from EpiGeEC comparing only user samples.(Please note that currently comparing
        large reference database with user samples is not supported).
        Finally, the pipeline executes four consecutive report steps to create the final report of the pipeline with quality control labels.

        This new pipeline can be used for pre-validation in order to assess the usability of a dataset in any given study, even
        in the absence of the original raw reads files. This is an advantage, for instance, in the case of
        human epigenomic datasets within IHEC, as signal tracks are made publicly available, while raw
        data files are stored in controlled access repositories.

        You can test this pipeline with ChIP-Seq samples from the IHEC portal :
        https://epigenomesportal.ca/ihec/grid.html?assembly=4&build=2018-10
    """

    def __init__(self, *args, protocol=None, **kwargs):
        if protocol is None:
            self._protocol = 'default'
        else:
            self._protocol = protocol
        super(EpiQC, self).__init__(*args, **kwargs)

    @property
    def output_dirs(self):
        dirs = {'bigwiginfo_output_directory': 'bigwiginfo',
                'chromimpute_output_directory': 'imputation',
                'bedgraph_converted_directory': 'bedgraph_data',
                'bedgraph_chr_converted_directory': 'bedgraph_chr_data',
                'chromimpute_converted_directory': 'converted',
                'chromimpute_distance_directory': 'distance',
                'chromimpute_traindata_directory': 'traindata',
                'chromimpute_predictor_directory': 'predictor',
                'chromimpute_apply': 'imputed',
                'chromimpute_eval': 'eval',
                'signal_to_noise_output_directory': 'signal_to_noise',
                'epigeec_output_directory': 'epigeec',
                'epigeec_hdf5': 'hdf5',
                'epigeec_filtered': 'filtered',
                'epigeec_output': 'output',
                'report_dir': 'report'
                }

        return dirs

    @property
    def inputinfo_file(self):
        inputinfo_filename = "inputinfofile.ChromImpute.txt"
        return inputinfo_filename

    @property
    def bigwiginfo_output(self):
        file_name = {'prefix': 'bigwiginfo',
                     'extension': '.txt'}
        return file_name

    @property
    def chipseq_bigwig(self):
        file_name = { 'tracks_dir' : 'tracks',
        'bigwig_dir': 'bigWig',
        'extension': '.bw'
        }
        return file_name

    @property
    def mark_type_conversion(self):
        dirs = {'N': 'narrow',
                'B': 'broad',
                'I': 'Input'
                }
        return dirs

    @property
    def readsets(self):
        if self._readsets is None:
            for readset in super().readsets:
                if not readset.mark_name:
                    _raise(SanitycheckError("Error: missing readset MarkName for " + readset.name))
                elif not readset.mark_type:
                    _raise(SanitycheckError("Error: missing readset MarkType for " + readset.name))
        return self._readsets

    @property
    def chromosome_file(self):
        file_name = os.path.join(self.output_dir, self.output_dirs['chromimpute_output_directory'],
                                      "_".join((global_config_parser.param('DEFAULT', 'scientific_name'), global_config_parser.param('DEFAULT',
                                                                                                         'assembly'),
                                                "chrom_sizes.txt")))
        return file_name

    @property
    def prefix_path(self):
        # Depending on the path of the readset, file this try to generate the path to singal tracks folder
        # Therefore, It is necessary to place the readset file in the ChIP-seq output directory
        #
        #gets the complete path to the readset directory and remove the readset name from the file path
        path = os.path.relpath(self.readsets_file.name, self.output_dir).replace(os.path.basename(self.readsets_file.name), '')

        return path

    def parseInputInfoFile(self, inputinfofile):
        """
            Parses a chromimpute input info file.
            Each element of the array samplesMarksFiles contains the name of the sample, its mark, and the name of the file in that order.
        """
        samplesMarksFiles = []
        for line in inputinfofile:
            sampleMarkFile = [line[0], line[1], line[2]]
            samplesMarksFiles.append(sampleMarkFile)

        return samplesMarksFiles

    def bigwiginfo(self):
        """
            Runs the tool bigWigInfo on bigwig files (

            Inspecting signal tracks to identify some obvious problems that
            could have an impact on the quality of the ChIP-Seq data is performed by UCSC-bigwiginfo
            https://bioconda-recipes-demo.readthedocs.io/en/docs/recipes/ucsc-bigwiginfo/README.html)
            bigWigInfo is capable of identifying obvious issues such as
            missing chromosomes and insufficient track coverage, which are usually symptoms of
            improperly generated tracks.

            If the user has specified bigwig files in the readset file under BIGIWIG column, they will be utilized by
            the tool. Otherwise, the user is required to process files using ChIp-Seq pipeline to generate
            bigwig files. Then paths for bigwig files are reconstructed based on the ChIP-Seq readset file
            and will be used subsequently.
            (Note that: the readset file should be in the same folder as the ChIp-Seq output)
        """
        jobs = []





        output_dir = self.output_dirs['bigwiginfo_output_directory']

        #obtain samples names for only histone marks (not inputs) from design
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    if readset.bigwig: # Check if the readset has a BIGWIG column != None
                        bigwig_file = readset.bigwig
                    else:                      # If not, we search for the path from a chipseq pipeline
                          # Find path to chipseq folder
                        bigwig_file = os.path.join(self.prefix_path, self.chipseq_bigwig['tracks_dir'], sample.name,
                                                   readset.mark_name, self.chipseq_bigwig['bigwig_dir'],
                                                   sample.name +"."+ readset.mark_name+ self.chipseq_bigwig['extension']) # Create path to bigwig file

                    output_file = os.path.join(output_dir,
                                          self.bigwiginfo_output['prefix'] + "_" + os.path.basename(bigwig_file) + self.bigwiginfo_output['extension'])
                    job = concat_jobs([
                            Job(command="mkdir -p " + output_dir),
                            bigwiginfo.bigWigInfo(bigwig_file, output_file)
                        ])
                    job.name = "bigwiginfo." + readset.name
                    job.samples = [readset.sample]
                    jobs.append(job)

        return jobs

    def bigwig_to_bedgraph(self):
        """
            ucsc-bigwigtobedgraph (https://bioconda-recipes-demo.readthedocs.io/en/docs/recipes/ucsc-bigwigtobedgraph/README.html)
            is used to convert bigwig files to bedgraph files (Used in ChromImpute subsequently).
        """
        jobs = []

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    if readset.bigwig:  # Check if the readset has a BIGWIG column != None
                        bigwig_file = readset.bigwig
                    else:  # If not, we search for the path from a chipseq pipeline
                          # Find path to chipseq folder
                        bigwig_file = os.path.join(self.prefix_path, self.chipseq_bigwig['tracks_dir'], sample.name,
                                                   readset.mark_name, self.chipseq_bigwig['bigwig_dir'],
                                                   sample.name + "." + readset.mark_name + self.chipseq_bigwig[
                                                       'extension'])  # Create path to bigwig file
                    output_bedgraph = os.path.join(self.output_dirs['bedgraph_converted_directory'],
                                                   sample.name + "_" + readset.mark_name + ".bedgraph")
                    output_bedgraph_gz = os.path.join(self.output_dirs['bedgraph_converted_directory'],
                                                      sample.name + "_" + readset.mark_name + ".bedgraph.gz")
                    # one job for each sample and each jobs has a job for create the folder if not exist
                    job = concat_jobs([
                        Job(command="mkdir -p " + self.output_dirs['bedgraph_converted_directory']),
                        bigwiginfo.bigWigToBedGraph(bigwig_file, output_bedgraph),
                        Job(output_files=[output_bedgraph_gz], command="gzip -f " + output_bedgraph)
                    ])
                    job.name = "bigwig_to_bedgraph." + sample.name + "_" + readset.mark_name
                    job.samples = [readset.sample]
                    jobs.append(job)

        return jobs

    def chromimpute_preprocess(self):
        """
            This step (mandatory) is performed to create chromimpute directories, chromosome sizes file, inputinfo file with IHEC
            and user samples and finally link the converted IHEC bedgraph files to user directory. In order to run
            the chromimpute, inputinfo and chromsizes file are required to be in the imputation directory.
            Please note: chromsizes and inputinfo files are created dynamically when the user runs the pipeline. It is
            not necessary to submit jobs to create them.


        """

        # for now there is no way we can start from bedgraph files without creating a folder called "bedgraph_data"
        # and put all the files there
        # there should be a way to define it in the readset file. or use a specific readset file for epiqc
        # for now pipeline does not support that
        jobs = []
        inputinfofile = os.path.join(self.output_dir, self.output_dirs['chromimpute_output_directory'], self.inputinfo_file)

        bedgraph_converted_files = []

        self.create_chr_sizes()
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    bedgraph_converted_files.append(
                        os.path.join(self.output_dirs['bedgraph_converted_directory'],
                                     sample.name + "_" + readset.mark_name + ".bedgraph.gz"))

        job_folder_create = Job(input_files=bedgraph_converted_files,
                                command="""\
mkdir -p \\
{output_dir}/{converteddir} \\
{output_dir}/{compute_global_dist} \\
{output_dir}/{generate_train_data} \\
{output_dir}/{train} \\
{output_dir}/{apply} \\
{output_dir}/{eval}""".format(
                                    output_dir=self.output_dirs['chromimpute_output_directory'],
                                    converteddir=self.output_dirs['chromimpute_converted_directory'],
                                    compute_global_dist=self.output_dirs['chromimpute_distance_directory'],
                                    generate_train_data=self.output_dirs['chromimpute_traindata_directory'],
                                    train=self.output_dirs['chromimpute_predictor_directory'],
                                    apply=self.output_dirs['chromimpute_apply'],
                                    eval=self.output_dirs['chromimpute_eval'])
                                )

        # load inputinfo file path from ini file. this file is stored in CVMFS and currently, file stored in user directory is not supported
        # because there is a prefix file path which called CVMFS directory
        # in future when need to add user inputinfo file. add another paramter in ini file and check it with a IF statement
        # to check whether user needs to add inputinfo file in epiqc.py
        # accordingly call SVMFS or user file

        #ihec_inputinfofile = "/home/pubudu/projects/rrg-bourqueg-ad/pubudu/epiqc_test/test_version2/inputinfofile.txt"
        ihec_inputinfofile = os.environ[global_config_parser.param('DEFAULT', 'mugqic_path')] + "/" + global_config_parser.param('chromimpute',
                                                                                                     'IHEC_inputinfo')

        #remove inputinfor file if exist
        #at this point imputation directory has already been created as chromosome file has generated first
        if os.path.exists(inputinfofile):
            os.remove(inputinfofile)
        #copy inputinfor file from CVMFS
        copyfile(ihec_inputinfofile, inputinfofile)

        #add user histone marks to inputinfo file
        # dynamically extend inputinfor file adding user samples (this is not a job, this step is executed
        # when generating the job script
        # copy histone mark, sample and file paths in user's samples into inputinfo file (avoid input histone files)
        with open(inputinfofile, "a") as inputinfo:
            for sample in self.samples:
                for readset in sample.readsets:
                    if readset.mark_type != "I":
                        inputinfo.write("{sample}\t{histone}\t{file}\n".format(
                            sample=sample.name,
                            histone=readset.mark_name,
                            file=sample.name + "_" + readset.mark_name + ".bedgraph.gz"))


        #train_user_data = config.param('DEFAULT', 'train_only_user_data')
        train_user_data = "F"
        chr_sizes_file = self.chromosome_file

        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_converted_directory'])

        inputinfofile = os.path.join(self.output_dir, self.output_dirs['chromimpute_output_directory'], self.inputinfo_file)

        #gather converted file paths in CVMFS as they needed to feed into the job as input files
        converted_simlinks = []
        if train_user_data == "F":
            with open(chr_sizes_file, "r") as chrominfofile:
                for line in chrominfofile:
                    chr_name = line.strip().split("\t")[0]
                    with open(inputinfofile, "r") as inputinfo:
                        for inputinfoline in inputinfo:
                            # if histone mark in inputinfo file present in readset file it only includes as an input fileof the convert global distance step

                            converted_simlinks.append(os.path.join(converteddir,
                                                                "%s_%s.wig.gz" %
                                                                (chr_name,
                                                                 inputinfoline.strip().split("\t")[2])))


        job = []
        #add chr_sizes_file and inputinfo file to list of input files
        converted_simlinks.extend([chr_sizes_file, inputinfofile])

        #job to create simlinks for converted signal tracks in IHEC data set stored in CVMFS
        #these will be used in ChromImpute
        job_create_simlinks = Job(
            output_files=converted_simlinks,
            command="""\
    if [ "$(ls -A {output_dir}/{user_converteddir})" ]; then
    rm {output_dir}/{user_converteddir}/*
    fi &&
    ln -s {ihec_converteddir}/* {output_dir}/{user_converteddir}/""".format(
            output_dir=self.output_dirs['chromimpute_output_directory'],
            user_converteddir=self.output_dirs['chromimpute_converted_directory'],
            ihec_converteddir=global_config_parser.param('chromimpute', 'IHEC_data') )
        )


        jobs.append(concat_jobs([job_folder_create, job_create_simlinks ],
                                name="chromimpute_preprocess"))

        return jobs

    def create_chr_sizes(self):

        output_dir = os.path.join(self.output_dir, self.output_dirs['chromimpute_output_directory'])

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        #get environment variable to generate the path of the chromosome file. different from usual ini file path.
        #check the ini file
        chr_sizes = os.environ[global_config_parser.param('DEFAULT', 'mugqic_path')] + "/" + global_config_parser.param('DEFAULT', 'chromosome_size')
        chrs = global_config_parser.param('chromimpute_preprocess', 'chromosomes')
        # get the chromosome from the ini file if one chr specified it gets the chromosome and split the string in the
        # ini file
        # otherwise get all chrs information from dictionary file and creates the chromosome length file
        if chrs == "All":
            genome_dict = os.path.expandvars(global_config_parser.param('DEFAULT', 'genome_dictionary', param_type='filepath'))
            chrs = genome.chr_names_conv(genome_dict)
        else:
            chrs = global_config_parser.param('chromimpute_preprocess', 'chromosomes').split(",")

        chr_sizes_file = self.chromosome_file

        if os.path.exists(chr_sizes_file):
            os.remove(chr_sizes_file)
        # Dynamically creates the chromosome size file using ini file parameters
        for chr in chrs:
            with open(chr_sizes, "r") as chr_sizes_genome:
                for chr_line in chr_sizes_genome:
                    if chr == chr_line.strip().split("\t")[0]:
                        chr_size = chr_line.strip().split("\t")[1]
                        chr_sizes_file_name = open(chr_sizes_file, "a")
                        chr_sizes_file_name.write("%s\t%i\n" % (chr, int(chr_size)))

    def chromimpute_convert(self):
        """
            This step is performed to convert each unique mark and sample combination signal tracks (bedgraph files)
            in the user dataset into binned signal resolution tracks.

            If you got index out of bound exception, check whether your reference genome version of the bedgraph file
            is similar to chromosome_sizes_file inside the imputation folder
        """
        jobs = []

        inputinfofile = os.path.join(self.output_dir, self.output_dirs['chromimpute_output_directory'], self.inputinfo_file)

        # check ini file whether user has requested specific chromosoems instead All chromosomes. Then change the input folder accordingly
        chrs = global_config_parser.param('chromimpute_preprocess', 'chromosomes')

        if chrs == "All":
            genome_dict = os.path.expandvars(global_config_parser.param('DEFAULT', 'genome_dictionary', param_type='filepath'))
            all_chrs = genome.chr_names_conv(genome_dict)
        else:
            all_chrs = global_config_parser.param('chromimpute_preprocess', 'chromosomes').split(",")

        input_dir = self.output_dirs['bedgraph_converted_directory']
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                  self.output_dirs['chromimpute_converted_directory'])
        chr_sizes_file = self.chromosome_file

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    output_files = []
                    for chr in all_chrs:
                        output_files.append(os.path.join(output_dir, chr + "_" + sample.name + "_" + readset.mark_name +".bedgraph.gz.wig.gz"))

                    input_file = os.path.join(input_dir,
                                     sample.name + "_" + readset.mark_name + ".bedgraph.gz")

                    job = chromimpute.convert(input_dir, input_file, output_dir, output_files, inputinfofile,
                                              readset.mark_name, sample.name, chr_sizes_file)
                    job.samples = [sample]
                    jobs.append(job)
        return jobs

    def chromimpute_compute_global_dist(self):
        """

            This steps is used to Compute the global distance based on correlation for each mark in each sample with
            the same mark in all other samples in inputinfo file. Creates a file for each mark in each sample containing
            a ranked list of the globally nearest samples.

        """
        jobs = []
        chr_sizes_file = self.chromosome_file

        inputinfofile = os.path.join(self.output_dir, self.output_dirs['chromimpute_output_directory'], self.inputinfo_file)
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_converted_directory'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                  self.output_dirs['chromimpute_distance_directory'])

        histone_marks = []
        # create a job for each histone maek in inputinfor file together with user's histone marks
        # get unique histone marks from the inputinfo file
        with open(inputinfofile, "r") as histone_marks_total:
            for line in histone_marks_total:
                histone_marks.append(line.strip().split("\t")[1])
        histone_marks = list(set(histone_marks))



        # this check all the converted files for converted step as input file
        # so it takes time
        # if you want to do this uncomment this and comment the code below this
        # for histone in histone_marks:
        #     with open(inputinfofile, "r") as inputinfo:
        #         for inputinfoline in inputinfo:
        #             #if histone mark in inputinfo file present in design file it only includes as an input fileof the convert global distance step
        #             if inputinfoline.split("\t")[1]==histone:
        #
        #                 with open(chr_sizes_file, "r") as chrominfofile:
        #                     for line in chrominfofile:
        #
        #                         input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'],
        #                                     self.output_dirs['chromimpute_converted_directory'], "%s_%s.wig.gz" %
        #                                             (line.strip().split("\t")[0], inputinfoline.strip().split("\t")[2])))
        #
        #                 output_files.append(
        #                     os.path.join(output_dir, "%s_%s.txt" % (inputinfoline.split("\t")[0], histone)))
        #
        #     job = chromimpute.compute_global_dist(input_files, output_dir, output_files, converteddir, inputinfofile,
        #                                           histone, chr_sizes_file)

        # this only get first line in the chromosome file(i.e chr1) as inputs to the converted step. if u want to get
        # all the chrms uncomment above comment and comment below code

        #don't need to add self.output_dir as the prefix here since they are jobs and job system can identify the
        # correct path
        for histone in histone_marks:
            input_files = []
            output_files = []
            with open(chr_sizes_file, "r") as chrominfofile:
                line = chrominfofile.readline()

                with open(inputinfofile, "r") as inputinfo:
                    for inputinfoline in inputinfo:

                        # if histone mark in inputinfo file present in design file it only includes as an input fileof the convert global distance step
                        if inputinfoline.split("\t")[1] == histone:

                            input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'],
                                                            self.output_dirs['chromimpute_converted_directory'],
                                                            "%s_%s.wig.gz" %
                                                            (line.strip().split("\t")[0],
                                                             inputinfoline.strip().split("\t")[2])))

                            output_files.append(
                                os.path.join(output_dir, "%s_%s.txt" % (inputinfoline.split("\t")[0], histone)))

            job = chromimpute.compute_global_dist(input_files, output_dir, output_files, converteddir, inputinfofile,
                                                  histone, chr_sizes_file)
            jobs.append(job)

        return jobs

    def chromimpute_generate_train_data(self):
        """
            This step is performed to generate a set of training data instances taking directory of converted data
            and global distances.
        """
        jobs = []

        # since the inputinfo file is called when generating job scripts, the prefix to identify the working
        # directory is necessary.

        inputinfofile = os.path.join(self.output_dir, self.output_dirs['chromimpute_output_directory'], self.inputinfo_file)
        temp2_inputinfofile_path = os.path.join(self.output_dir, self.output_dirs['chromimpute_output_directory'],
                                                "temp2_" + self.inputinfo_file)
        temp_inputinfofile_path = os.path.join(self.output_dir, self.output_dirs['chromimpute_output_directory'],
                                               "temp_" + self.inputinfo_file)
        distancedir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                   self.output_dirs['chromimpute_distance_directory'])
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_converted_directory'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                  self.output_dirs['chromimpute_traindata_directory'])
        chr_sizes_file = self.chromosome_file
        #get all the unique histone marks in user's readset file

        histone_marks =[]
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    histone_marks.append(readset.mark_name)

        histone_marks = set(histone_marks)

        # train_data_path = config.param('chromimpute_generate_train_data', 'pre_trained_data_path') #currently not supported
        train_data_path = ""
        # create temp inputinfo file with interneal indexes used in traindata step
        # this file will be used to get the integer index for list of output files

        # job_temp_inputinfo = chromimpute.temp_inputinfo( inputinfofile, temp_inputinfofile)
        # jobs.append(job_temp_inputinfo)

        # train_user_data = config.param('DEFAULT', 'train_only_user_data') #currently not supported
        train_user_data = "F"
        # This file is generated by python dynamically while runnning GenPipes.
        # Note that, the user does not need to submit the job script file to the server to create this file.
        # Since GenPipes amend text in the file, removing this file is essential everytime when running GenPipes.
        # We cannot use tmpdir or bash to remove this file and need it to do by GenPipes. So adding os.remove is intentional
        if train_user_data == "F":
            if os.path.exists(temp2_inputinfofile_path):
                os.remove(temp2_inputinfofile_path)

            if os.path.exists(temp_inputinfofile_path):
                os.remove(temp_inputinfofile_path)

            temp_inputinfofile = open(temp2_inputinfofile_path, "a")
            lines = open(inputinfofile, 'r').readlines()

            for line in sorted(lines, key=lambda line: line.split()[0]):
                temp_inputinfofile.write(line)
            temp_inputinfofile.close()
            temp_inputinfofile = open(temp_inputinfofile_path, "a")
            with open(temp2_inputinfofile_path, "r") as inputinfo:
                i = 1
                index_i = 0
                for inputinfoline in inputinfo:
                    if (i == 1):
                        sample = inputinfoline.split("\t")[0]
                        index_i = i - 1
                        temp_inputinfofile.write("%s\t%i%s" % (inputinfoline.strip(), index_i, os.linesep))
                        i += 1
                    else:
                        if sample == inputinfoline.split("\t")[0]:
                            temp_inputinfofile.write("%s\t%i%s" % (inputinfoline.strip(), index_i, os.linesep))
                            i += 1
                        else:
                            sample = inputinfoline.split("\t")[0]
                            index_i += 1
                            temp_inputinfofile.write("%s\t%i%s" % (inputinfoline.strip(), index_i, os.linesep))
                            i += 1
            os.remove(temp2_inputinfofile_path)
            temp_inputinfofile.close()
            if train_data_path == "":

                #        jobs.append(chromimpute.generate_train_data(input_files, output_dir, output_files, converteddir,
                #        distancedir, inputinfofile, contrast.real_name))

                for histone in histone_marks:
                    with open(chr_sizes_file, "r") as chrominfofile:
                        for line in chrominfofile:
                            chr_name = line.strip().split("\t")[0]
                            input_files = []
                            output_files = []
                            with open(temp_inputinfofile_path, "r") as inputinfo:

                                for inputinfoline in inputinfo:

                                    # if histone mark in inputinfo file present in readset file it only includes as an input fileof the convert global distance step
                                    if inputinfoline.split("\t")[1] == histone:
                                        input_files.append(os.path.join(distancedir, "%s_%s.txt" % (
                                            inputinfoline.split("\t")[0], histone)))
                                        # to make this job step independent of global distance step we need to correctly
                                        # specify all the dependency input files
                                        with open(inputinfofile, "r") as inputinfo2:
                                            for inpputinfo2line in inputinfo2:
                                                input_files.append(os.path.join(distancedir, "%s_%s.txt" % (
                                                inpputinfo2line.split("\t")[0], inpputinfo2line.split("\t")[1])))
                                        # distance files generated for user histone marks
                                        input_files.append(os.path.join(converteddir,
                                                                        "%s_%s.wig.gz" %
                                                                        (chr_name,
                                                                         inputinfoline.strip().split("\t")[2])))

                                        output_files.append(
                                            os.path.join(output_dir, "%s_traindata_%s_%i_0.txt.gz" % (
                                                chr_name, histone, int(inputinfoline.strip().split("\t")[3]))))
                                        output_files.append(
                                            os.path.join(output_dir, "attributes_%s_%i_0.txt.gz" % (
                                                histone, int(inputinfoline.strip().split("\t")[3]))))


                    # input_files.append(temp_inputinfofile)
                    job = chromimpute.generate_train_data(input_files, output_dir, output_files, converteddir,
                                                          distancedir,
                                                          inputinfofile, histone, chr_sizes_file, chr_name)
                    jobs.append(job)

            else:
                job_create_simlinks = Job(command="""\
                    if [ "$(ls -A {output_dir}/{usertraindatadir})" ]; then
                    rm {output_dir}/{usertraindatadir}/*
                    fi &&
                    ln -s {ihec_traindatadir}/* {output_dir}/{usertraindatadir}/""".format(
                    output_dir=self.output_dirs['chromimpute_output_directory'],
                    usertraindatadir=self.output_dirs['chromimpute_traindata_directory'],
                    ihec_traindatadir=train_data_path)
                )
                jobs.append(
                    concat_jobs([job_create_simlinks], name="chromimpute_generate_train_data.use_pre_trained_data"))
        else:
            log.info("Currently, only training user data is not supported. Skipping ....")
            # todo
            # mula idanma hadan enna one meka
            # user data witharak train karanna

        # chroms = []
        # with open(os.path.join(os.environ[config.param('chromimpute', 'chromosome_size').split("/")[0].replace("$", "")], config.param('chromimpute', 'chromosome_size').replace(config.param('chromimpute', 'chromosome_size').split("/")[0]+"/", "")), "r") as chrominfofile:#with open(config.param('chromimpute', 'chrominfofile')) as chrominfofile:             os.environ[config.param('chromimpute', 'chromosome_size').split("/")[0].replace("$", "")]
        #     for line in chrominfofile:
        #         chroms.append(line.split("\t")[0])

        # for contrast in self.contrasts:
        #     if contrast.treatments:
        #         input_files = []
        #         output_files = []
        #         indice = 0
        #         for sample in contrast.treatments:
        #             with open(os.path.join(os.environ[config.param('chromimpute', 'chromosome_size').split("/")[0].replace("$", "")], config.param('chromimpute', 'chromosome_size').replace(config.param('chromimpute', 'chromosome_size').split("/")[0]+"/", "")), "r") as chrominfofile:
        #                 for line in chrominfofile:
        #                     input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converted_directory'], "%s_%s.bedgraph.gz.wig.gz" % (line.split("\t")[0], sample.name)))
        #             input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_distance_directory'], "%s_%s.txt" % (sample.name, contrast.real_name)))
        #             output_files.append(os.path.join(output_dir, "attributes_%s_%i_0.txt.gz" % (contrast.real_name, indice)))
        #             output_files.append(os.path.join(output_dir, "traindata_%s_%i_0.txt.gz" % (contrast.real_name, indice)))
        #             indice += 1
        #         # for chrom in chroms:
        #         jobs.append(chromimpute.generate_train_data(input_files, output_dir, output_files, converteddir, distancedir, inputinfofile, contrast.real_name))

        return jobs

    def chromimpute_train(self):
        """
            This step is used to train regression trees based on the feature data produced by GenerateTrainData.

        """
        jobs = []

        temp_inputinfofile = os.path.join(self.output_dir, self.output_dirs['chromimpute_output_directory'],
                                          "temp_" + self.inputinfo_file)
        inputinfofile = os.path.join(self.output_dir, self.output_dirs['chromimpute_output_directory'],
                                     self.inputinfo_file)
        traindatadir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_traindata_directory'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                  self.output_dirs['chromimpute_predictor_directory'])

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    input_files = []
                    output_files = []
                    with open(temp_inputinfofile, "r") as inputinfo:
                        for inputinfoline in inputinfo:
                            inputinfo_histone = inputinfoline.strip().split("\t")[1]
                            inputinfo_sample = inputinfoline.strip().split("\t")[0]
                            if readset.mark_name == inputinfo_histone:
                                if sample.name != inputinfo_sample:
                                    output_files.append(os.path.join(output_dir, "useattributes_%s_%s_%i_0.txt.gz" % (
                                        sample.name, readset.mark_name, int(inputinfoline.strip().split("\t")[3]))))

                                input_files.append(
                                    os.path.join(traindatadir, "attributes_%s_%i_0.txt.gz" % (
                                        readset.mark_name, int(inputinfoline.strip().split("\t")[3]))))

                    jobs.append(
                    chromimpute.train(input_files, output_dir, output_files, traindatadir, inputinfofile, sample.name,
                                      readset.mark_name))

        return jobs

    def chromimpute_apply(self):
        """
            This step is used to apply the predictors generated in the Train command to generate the imputed data.
            A job is created for each sample and mark given in the dataset.
        """
        jobs = []

        temp_inputinfofile = os.path.join(self.output_dir, self.output_dirs['chromimpute_output_directory'],
                                          "temp_" + self.inputinfo_file)
        inputinfofile = os.path.join(self.output_dir, self.output_dirs['chromimpute_output_directory'],
                                     self.inputinfo_file)
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_converted_directory'])
        distancedir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                   self.output_dirs['chromimpute_distance_directory'])
        predictordir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_predictor_directory'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                  self.output_dirs['chromimpute_apply'])
        chr_sizes_file = self.chromosome_file

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    input_files = []
                    output_files = []
                    with open(chr_sizes_file, "r") as chrominfofile:
                        for line in chrominfofile:
                            chr_name = line.strip().split("\t")[0]
                            with open(temp_inputinfofile, "r") as inputinfo:
                                for inputinfoline in inputinfo:
                                    inputinfo_histone = inputinfoline.strip().split("\t")[1]
                                    inputinfo_sample = inputinfoline.strip().split("\t")[0]
                                    if readset.mark_name == inputinfo_histone:
                                        if sample.name != inputinfo_sample:
                                            input_files.append(
                                                os.path.join(predictordir, "useattributes_%s_%s_%i_0.txt.gz" % (
                                                    sample.name, readset.mark_name,
                                                    int(inputinfoline.strip().split("\t")[3]))))

                                # output_files.append(os.path.join(output_dir, "classifier_%s_%s_%i_0.txt.gz" % (
                                #   sample.name, contrast.real_name, int(inputinfoline.strip().split("\t")[3]))))
                            output_files.append(
                                os.path.join(output_dir, "%s_impute_%s_%s.wig.gz" % (chr_name, sample.name,
                                                                                     readset.mark_name)))

                            jobs.append(
                                chromimpute.apply(input_files, output_dir, converteddir, distancedir, predictordir,
                                                  inputinfofile, sample.name, readset.mark_name, chr_name,
                                                  chr_sizes_file, output_files))

        return jobs

    def chromimpute_eval(self):
        """
            This step is used to compare the agreement between an observed and imputed data set.
            A job is created for every sample-mark given in the dataset.
        """
        jobs = []

        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_converted_directory'])
        imputeddir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                  self.output_dirs['chromimpute_apply'])
        chr_sizes_file = self.chromosome_file

        percent1 = global_config_parser.param('chromimpute_eval', 'percent1')
        percent2 = global_config_parser.param('chromimpute_eval', 'percent2')
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs[
            'chromimpute_eval'])

        for sample in self.samples:
            for readset in sample.readsets:
                input_files = []
                if readset.mark_type != "I":
                    with open(chr_sizes_file, "r") as chrominfofile:
                        line = chrominfofile.readline()
                        chr_name = line.strip().split("\t")[0]

                        input_files.append(os.path.join(imputeddir, "%s_impute_%s_%s.wig.gz" % (
                        chr_name, sample.name, readset.mark_name)))

                        input_files.append(os.path.join(converteddir, "%s_%s_%s.bedgraph.gz.wig.gz" % (
                        chr_name, sample.name, readset.mark_name)))

                        imputed_file = "impute_%s_%s.wig.gz" % (
                        sample.name, readset.mark_name)

                        converted_file = "%s_%s.bedgraph.gz.wig.gz" % (
                        sample.name, readset.mark_name)

                    output_file = os.path.join(output_dir, "eval_%s_%s.tsv" % (sample.name,
                                                                               readset.mark_name))

                    jobs.append(chromimpute.eval(input_files, imputed_file, converted_file, output_file, converteddir,
                                                  imputeddir, percent1, percent2, chr_sizes_file, sample.name,
                                                  readset.mark_name))

        return jobs

    def chromimpute(self):

        """
        Runs the steps (bigwig_to_bedgraph, chromimpute_preprocess, Convert, ComputeGlobalDist,
        GenerateTrainData, Train apply and eval) of the ChromImpute on the bigwig files given to the pipeline.
            All the output are stored in the imputation directory.

        """
        jobs = []
        jobs.extend(self.bigwig_to_bedgraph())
        jobs.extend(self.chromimpute_preprocess())
        jobs.extend(self.chromimpute_convert())
        jobs.extend(self.chromimpute_compute_global_dist())
        jobs.extend(self.chromimpute_generate_train_data())
        jobs.extend(self.chromimpute_train())
        jobs.extend(self.chromimpute_apply())
        jobs.extend(self.chromimpute_eval())
        return jobs

    def signal_to_noise(self):
        """
            Binned signal resolution tracks from chromipute convert step is used to determined the
            percentage of the whole file signal that was located in the top 5% and 10% of the bins
            The resulted output is a tsv file
        """
        jobs = []
        output_dir = self.output_dirs['signal_to_noise_output_directory']
        chr_sizes_file = self.chromosome_file

        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_converted_directory'])

        self.create_chr_sizes()

        for sample in self.samples:
            for readset in sample.readsets:
                input_files = []
                if readset.mark_type != "I":
                    with open(chr_sizes_file, "r") as chrominfofile:
                        for line in chrominfofile:
                            chr_name = line.strip().split("\t")[0]
                            # If not, we search for the path from a chipseq pipeline
                            converted_bedgraph_file = os.path.join(converteddir,
                                                                   "%s_%s_%s.bedgraph.gz.wig.gz" % (chr_name, sample.name, readset.mark_name))

                            output_file = os.path.join(output_dir,
                                                       "%s_%s_%s.tsv" % (chr_name, sample.name, readset.mark_name))
                            signal_noise_folder = Job(
                                command="mkdir -p {output_dir}".format(
                                    output_dir=output_dir
                                ))
                            signal_noise_job = Job(
                                [converted_bedgraph_file, chr_sizes_file],
                                [output_file],
                                [['signal_noise', 'module_python'],
                                 ['signal_noise', 'module_mugqic_tools']],
                                command="""\
python $PYTHON_TOOLS/signal_noise.py \\                               
-i {input_file} \\
-p1 {percent1} \\
-p2 {percent2} \\
-o {output_dir}""".format(
                                    input_file=converted_bedgraph_file,
                                    percent1=global_config_parser.param('signal_noise', 'percent1'),
                                    percent2=global_config_parser.param('signal_noise', 'percent2'),
                                    output_dir=output_file
                                ))

                            jobs.append(concat_jobs([signal_noise_folder, signal_noise_job], name='signal_noise.' + "%s_%s_%s.tsv" % (chr_name, sample.name, readset.mark_name)))

        return jobs

    def epigeec_tohdf5(self):

        #create hdf5 files
        jobs = []

        input_dir = 'epigeec'
        output_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_hdf5'])
        job_epigeec_hdf5 =[]
        hdf5_folder = Job(
            command="mkdir -p {output_dir}".format(
                output_dir=output_dir
            ))
        for sample in self.samples:
            for readset in sample.readsets:
                input_files = []
                if readset.mark_type != "I":
                    if readset.bigwig:  # Check if the readset has a BIGWIG column != None
                        bigwig_file = readset.bigwig
                    else:  # If not, we search for the path from a chipseq pipeline
                          # Find path to chipseq folder
                        bigwig_file = os.path.join(self.prefix_path, "tracks", sample.name, readset.mark_name, "bigWig",
                                                   sample.name + "." + readset.mark_name + ".bw")  # Create path to bigwig file

                    job_epigeec_hdf5.append(epigeec.tohdf5(output_dir, bigwig_file))
                    jobs_epigeec_hdf5 = concat_jobs(job_epigeec_hdf5)
        jobs.append(concat_jobs([hdf5_folder, jobs_epigeec_hdf5],name = "epigeec_tohdf5"))

        return jobs

    def epigeec_filter(self, hdf5_files):

        # filter regions from hdf5 files
        jobs = []
        filter_job = []
        input_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_hdf5'])
        output_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_filtered'])
        filter_folder_job = Job(
            command="mkdir -p {output_dir}".format(
                output_dir=output_dir
            ))



        for hdf5_file in hdf5_files:
            path_to_hdf5_file = os.path.join(input_dir, hdf5_file)
            filter_job.append(epigeec.filter(output_dir, path_to_hdf5_file))
            filter_jobs = concat_jobs(filter_job)

        jobs.append(concat_jobs([filter_folder_job, filter_jobs],name="epigeec_filter"))

        return jobs

    def epigeec_correlate(self, skip_filter_step, hdf5_files):

        #create the correlate matrix file
        jobs=[]
        job=[]
        hdf5_file_list_name = os.path.join(self.output_dirs['epigeec_output_directory'],
                                                 "hdf5_list.txt")

        input_files =[]
        output_dir = os.path.join(self.output_dirs['epigeec_output_directory'],
                                                 self.output_dirs['epigeec_output'])
        output_folder_job = Job(
            command="mkdir -p {output_dir}".format(
                output_dir=output_dir
            ))

        job_create_hdf5_list_file = Job(
            output_files=[hdf5_file_list_name],
            command="""\
        if test -f "{hdf5_file_list_name}"; then
        rm {hdf5_file_list_name} &&
        touch {hdf5_file_list_name}
        else
            touch {hdf5_file_list_name}
        fi""".format(
            hdf5_file_list_name=hdf5_file_list_name)
        )

        for hdf5_file in hdf5_files:
            if skip_filter_step:
                path_to_hdf5_file = os.path.join(self.output_dirs['epigeec_output_directory'],
                                                 self.output_dirs['epigeec_hdf5'], hdf5_file)
                input_files.append(path_to_hdf5_file)
            else:
                path_to_hdf5_file = os.path.join(self.output_dirs['epigeec_output_directory'],
                                                 self.output_dirs['epigeec_filtered'], hdf5_file)
                input_files.append(path_to_hdf5_file)

            job_hdf5_list_file = epigeec.generate_hdf5_list(hdf5_file_list_name, path_to_hdf5_file)
            job.append(job_hdf5_list_file)
            jobs_hdf5_list_file = concat_jobs(job)


        job_matrix = epigeec.correlate(input_files, output_dir, hdf5_file_list_name)
        jobs.append(concat_jobs([output_folder_job, job_create_hdf5_list_file, jobs_hdf5_list_file, job_matrix], name="epigeec_correlate"))

        return jobs

    def epigeec(self):

        """
            Runs the epigeec pipeline (https://bitbucket.org/labjacquespe/epigeec/src/master/)
            on the bigwig files given to the pipeline
            Epigeec pipeline is consisted of 3 sub-steps
            1. Bigwig files are first converted to the hdf5 format
            2. Filter or select provided regions as a bed file (include or exclude) [optional]. The user can specify
            the options and the brf file path in the ini file. Otherwise, this step will be skipped
            3. Finally the correlation matrix is computed

        """
        jobs = []
        #check wjether filter files are specified in the ini, if so get the paths
        #paths should be specified as absolute paths
        #if not filter step will be skipped
        skip_filter_step = global_config_parser.param('epigeec', 'select') == '' and global_config_parser.param('epigeec', 'exclude') == ''

        hdf5_files = []

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    if readset.bigwig:  # Check if the readset has a BIGWIG column != None
                        bigwig_file = readset.bigwig
                    else:  # If not, we search for the path from a chipseq pipeline
                        prefix_path = "/".join(self.readsets_file.name.split("/")[:-1])  # Find path to chipseq folder
                        bigwig_file = os.path.join(prefix_path, "tracks", sample.name, readset.mark_name, "bigWig",
                                                   sample.name + "." + readset.mark_name + ".bw")  # Create path to bigwig file

                    hdf5_files.append(os.path.basename(bigwig_file) + ".hdf5")

        jobs.extend(self.epigeec_tohdf5())

        if not skip_filter_step:  # We skip this step if there are no filter files

            jobs.extend(self.epigeec_filter(hdf5_files))


        jobs.extend(self.epigeec_correlate(skip_filter_step, hdf5_files))

        return jobs

    def bigwig_info_report(self):

        """
        This step is performed to generate report on bigwiginfo result
        """
        jobs = []
        job = []
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    if readset.bigwig:  # Check if the readset has a BIGWIG column != None
                        bigwig_file = readset.bigwig

                    else:  # If not, we search for the path from a chipseq pipeline
                        prefix_path = "/".join(
                            self.readsets_file.name.split("/")[:-1])  # Find path to chipseq folder
                        bigwig_file = os.path.join(prefix_path, self.chipseq_bigwig['tracks_dir'], sample.name,
                                                   readset.mark_name, self.chipseq_bigwig['bigwig_dir'],
                                                   sample.name + "." + readset.mark_name + self.chipseq_bigwig[
                                                       'extension'])  # Create path to bigwig file

                    bigwiginfo_file = os.path.join(self.output_dirs['bigwiginfo_output_directory'],
                                                   self.bigwiginfo_output['prefix'] + "_" + os.path.basename(
                                                       bigwig_file) +
                                                   self.bigwiginfo_output['extension'])
                    report_file = os.path.join(self.output_dirs['report_dir'], sample.name,
                                                               readset.mark_name,
                                               "BigWigInfo_report_" + sample.name + "_" + readset.mark_name + ".txt")

                    bigwig_job = concat_jobs([
                        Job(command="mkdir -p " + os.path.join(self.output_dirs['report_dir'], sample.name,
                                                               readset.mark_name)),
                        epiqc_reports.bigwiginfo_report(bigwiginfo_file, report_file)
                    ])

                    job.append(bigwig_job)
                    bigwiginfo_report = concat_jobs(job)

        job = concat_jobs([  bigwiginfo_report])
        job.samples = self.samples
        job.name = "epiqc_report.bigwig_info"
        jobs.append(job)
        return jobs

    def signal_to_noise_report(self):
        """
        This step is performed to generate report on signal_to_noise result
        """
        jobs = []
        chr_job = []

        chr_sizes_file = self.chromosome_file
        report_folder_job = Job(
            command="mkdir -p {output_dir}".format(
                output_dir=self.output_dirs['report_dir']
            ))

        signal_noise_output_dir = self.output_dirs['signal_to_noise_output_directory']

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":

                    report_file = os.path.join(self.output_dirs['report_dir'], sample.name, readset.mark_name,
                                               "SignalToNoise_report_" + sample.name + "_" + readset.mark_name + ".txt")

                    with open(chr_sizes_file, "r") as chrominfofile:
                        for line in chrominfofile:
                            chr_name = line.strip().split("\t")[0]
                            signal_noise_file = os.path.join(signal_noise_output_dir,
                                                             "%s_%s_%s.tsv" % (
                                                             chr_name, sample.name, readset.mark_name))
                            #run through all the chrs since a signal_to_noise file has created for each chr

                            signalnoise_job = concat_jobs([
                                Job(command="mkdir -p " + os.path.join(self.output_dirs['report_dir'], sample.name,
                                                                       readset.mark_name)),
                                epiqc_reports.signal_to_noise_report(signal_noise_file, report_file, chr_name)
                            ])

                            chr_job.append(signalnoise_job)

                            signal_to_noise_report = concat_jobs(chr_job)

        job = concat_jobs([ signal_to_noise_report])
        job.samples = self.samples
        job.name = "epiqc_report.signal_to_noise"
        jobs.append(job)
        return jobs

    def chromimpute_report(self):
        """
        This step is performed to generate a report comparing ChromImpute imputed signal
track and input signal track (in bedgraph format).
        """
        jobs = []
        impute_job = []


        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":

                    report_file = os.path.join(self.output_dirs['report_dir'], sample.name, readset.mark_name,
                                               "ChromImpute_report_" + sample.name + "_" + readset.mark_name + ".txt")
                    eval_file = "eval_" + sample.name + "_" + readset.mark_name + ".tsv"
                    eval_file = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                             self.output_dirs['chromimpute_eval'], eval_file)

                    chroimpute_job = concat_jobs([
                        Job(command="mkdir -p " + os.path.join(self.output_dirs['report_dir'], sample.name,
                                                               readset.mark_name)),
                        epiqc_reports.chromimpute_report(eval_file, report_file)
                    ])

                    impute_job.append(chroimpute_job)

                    chromimpute_report = concat_jobs(impute_job)
        job = concat_jobs([chromimpute_report])
        job.samples = self.samples
        job.name = "epiqc_report.chromimpute"
        jobs.append(job)
        return jobs

    def epigeec_report(self):
        """
        This step is performed to generate a heatmap from EpiGeEC results
        """
        jobs = []

        report_folder_job = Job(
            command="mkdir -p {output_dir}".format(
                output_dir=self.output_dirs['report_dir']
            ))

        input_dir = os.path.join(self.output_dirs['epigeec_output_directory'],
                                 self.output_dirs['epigeec_output'])
        report_dir = os.path.join(self.output_dirs['report_dir'])
        output_heatmap_file = os.path.join(report_dir, "correlation_matrix.png")
        input_matrix_file = os.path.join(input_dir, "correlation_matrix.tsv")
        epigeec_report = epiqc_reports.epigeec_report(input_matrix_file, output_heatmap_file, report_dir)
        job = concat_jobs([report_folder_job, epigeec_report])
        job.samples = self.samples
        job.name = "epiqc_report.epigeec"
        jobs.append(job)
        return jobs


    def epiqc_report(self):

        jobs = []
        chr_job = []

        chr_sizes_file = self.chromosome_file

        signal_noise_output_dir = self.output_dirs['signal_to_noise_output_directory']

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    if readset.bigwig:  # Check if the readset has a BIGWIG column != None
                        bigwig_file = readset.bigwig

                    else:  # If not, we search for the path from a chipseq pipeline
                        prefix_path = "/".join(
                            self.readsets_file.name.split("/")[:-1])  # Find path to chipseq folder
                        bigwig_file = os.path.join(prefix_path, self.chipseq_bigwig['tracks_dir'], sample.name,
                                                   readset.mark_name, self.chipseq_bigwig['bigwig_dir'],
                                                   sample.name + "." + readset.mark_name + self.chipseq_bigwig[
                                                       'extension'])  # Create path to bigwig file

                    bigwiginfo_file = os.path.join(self.output_dirs['bigwiginfo_output_directory'],
                                                   self.bigwiginfo_output['prefix'] + "_" + os.path.basename(
                                                       bigwig_file) +
                                                   self.bigwiginfo_output['extension'])
                    report_file = os.path.join(self.output_dirs['report_dir'],
                                               "report_" + sample.name + "_" + readset.mark_name + ".txt")
                    report_folder_job = Job(
                        command="mkdir -p {output_dir}".format(
                            output_dir=self.output_dirs['report_dir']
                        ))

                    bigwiginfo_report = epiqc_reports.bigwiginfo_report(bigwiginfo_file,report_file)

                    with open(chr_sizes_file, "r") as chrominfofile:
                        for line in chrominfofile:
                            chr_name = line.strip().split("\t")[0]
                            signal_noise_file = os.path.join(signal_noise_output_dir,
                                                             "%s_%s_%s.tsv" % (
                                                             chr_name, sample.name, readset.mark_name))
                            #run through all the chrs since a signal_to_noise file has created for each chr
                            chr_job.append(epiqc_reports.signal_to_noise_report(signal_noise_file, report_file, chr_name))
                            signal_to_noise_report = concat_jobs(chr_job)

                    eval_file = "eval_" + sample.name + "_" + readset.mark_name + ".tsv"
                    eval_file = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                             self.output_dirs['chromimpute_eval'], eval_file)
                    chromimpute_report = epiqc_reports.chromimpute_report(eval_file, report_file)

                    job = concat_jobs(
                        [report_folder_job, bigwiginfo_report, signal_to_noise_report, chromimpute_report])
                    job.name = "epiqc_report." + readset.name + "_" + readset.mark_name
                    job.samples = self.samples
                    jobs.append(job)

        input_dir = os.path.join(self.output_dirs['epigeec_output_directory'],
                                 self.output_dirs['epigeec_output'])
        report_dir = os.path.join(self.output_dirs['report_dir'])
        output_heatmap_file = os.path.join(report_dir, "correlation_matrix.png")
        input_matrix_file = os.path.join(input_dir, "correlation_matrix.tsv")
        epigeec_report = epiqc_reports.epigeec_report(input_matrix_file, output_heatmap_file, report_dir)
        jobs.append(epigeec_report)
        return jobs

    def epiqc_final_report(self):

        """
               Creates a report file for each bigwig file.

               Alert levels :
                   High Level Alert:
                       Chromosome count is under 23
                   Medium Level Alert:
                       Whole genome bases covered under 25,000,000
                       GeEC average correlation score under 50% within the same consortium tracks
                       Signal in top 10% bins below 30%
                       ChromImpute OBSERVED_1.0_IMPUTE_5.0 below 30%
                   Low Level Alert:
                       Whole genome bases covered under 75,000,000
                       Signal in top 5% bins below 20%
                       ChromImpute BOTH_1.0 below 20%
           """
        jobs = []
        final_report = os.path.join(self.output_dirs['report_dir'],
                                    "epiqc_report" + ".tsv")
        report_jobs =[]

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    bigwiginfo_report = os.path.join(self.output_dirs['report_dir'], sample.name, readset.mark_name,
                                                      "BigWigInfo_report_" + sample.name + "_" + readset.mark_name + ".txt")
                    chromimpute_report = os.path.join(self.output_dirs['report_dir'], sample.name, readset.mark_name,
                                               "ChromImpute_report_" + sample.name + "_" + readset.mark_name + ".txt")
                    signalnoise_report = os.path.join(self.output_dirs['report_dir'], sample.name, readset.mark_name,
                                 "SignalToNoise_report_" + sample.name + "_" + readset.mark_name + ".txt")
                    report_job = concat_jobs([
                        epiqc_reports.final_report(bigwiginfo_report, chromimpute_report, signalnoise_report, final_report, sample.name, readset.mark_name)
                    ])
                    report_jobs.append(report_job)

                    epiqc_report = concat_jobs(report_jobs)


        report_file = Job(
            output_files=[final_report],
            command="""\
if test -f "{final_report}"; then
rm {final_report} && touch {final_report}
else
touch {final_report} &&
echo -e "Sample_name\tHistone_mark\tDecision" >> {final_report}
fi""".format(
                final_report=final_report)
        )
        job = concat_jobs([report_file, epiqc_report])
        job.samples = self.samples
        job.name = "epiqc_report.final_report"
        jobs.append(job)

        return jobs

    @property
    def step_list(self):
        return self.protocols()[self._protocol]

    def protocols(self):
        return { 'default': [
            self.bigwiginfo,
            self.chromimpute,
            self.signal_to_noise,
            self.epigeec,
            self.bigwig_info_report,
            self.chromimpute_report,
            self.signal_to_noise_report,
            self.epigeec_report,
            self.epiqc_final_report
            ] }

def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    # Check if Genpipes must be ran inside a container
    utils.container_wrapper_argparse(__file__, argv)
    # Build help
    epilog = EpiQC.process_help(argv)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve', epilog=epilog)

    # populate the parser
    parser = EpiQC.argparser(parser)

    parsed_args = parser.parse_args(argv)

    sanity_check = parsed_args.sanity_check
    loglevel = parsed_args.log
    utils.set_logger(loglevel, sanity_check=sanity_check)

    # Pipeline config
    config_files = parsed_args.config

    # Common Pipeline options
    genpipes_file = parsed_args.genpipes_file
    container = parsed_args.container
    clean = parsed_args.clean
    report = parsed_args.report
    no_json = parsed_args.no_json
    force = parsed_args.force
    job_scheduler = parsed_args.job_scheduler
    output_dir = parsed_args.output_dir
    steps = parsed_args.steps
    readset_file = parsed_args.readsets_file
    design_file = parsed_args.design_file


    pipeline = EpiQC(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file,
                         clean=clean, report=report, force=force, job_scheduler=job_scheduler, output_dir=output_dir,
                         design_file=design_file, no_json=no_json, container=container)

    pipeline.submit_jobs()

if __name__ == '__main__':
    main()
