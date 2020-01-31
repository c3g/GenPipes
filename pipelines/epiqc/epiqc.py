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
import os
import sys
import csv
from operator import attrgetter

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *

from pipelines import common

from bfx import bigwiginfo
from bfx import chromimpute
from bfx import wigSignalNoise
from bfx import epigeec
from bfx import epiqc_report
from bfx.readset import *

from pipelines.chipseq import chipseq

log = logging.getLogger(__name__)

class EpiQC(common.Illumina):
    """
        EpiQC Pipeline
        ==============

        EpiQC is a quality control pipeline for signal files (bigwig) generated from ChIP-Seq. The pipeline does a series of calculations on
        these files to determine wether or not they can be considered good.Four metrics are computed from a single bigwig file.
        With BigWigInfo it is determined if there are missing chromosomes, if the chromosome count is lower than 23 it raises a high level alert.
        ChromImpute imputes signal tracks for the first chromosome and using these imputed files EpiQC computes 2 other metrics.
        And finally the pipeline creates a heatmap from the correlation matrix obtained with EpiGeEC.

        You can test this pipeline with ChIP-Seq samples from the IHEC portal :
        https://epigenomesportal.ca/ihec/grid.html?assembly=4&build=2018-10
    """
    def __init__(self, protocol=None):
        self._protocol = protocol
        # Add pipeline specific arguments
        self.argparser.add_argument("-d", "--design", help="design file", type=file)
        super(EpiQC, self).__init__(protocol)

    @property
    def output_dirs(self):
        dirs = {'bigwiginfo_output_directory': 'bigwiginfo',
                'chromimpute_output_directory': 'imputation',
                'bedgraph_converted_directory' : 'bedgraph_data',
                'chromimpute_converted_directory' : 'converted',
                'chromimpute_distance_directory' : 'distance',
                'chromimpute_traindata_directory' : 'traindata',
                'chromimpute_predictor_directory' : 'predictor',
                'chromimpute_apply' : 'imputed',
                'chromimpute_eval' : 'eval',
                'signal_to_noise_output_directory': 'signal_to_noise',
                'epigeec_output_directory': 'epigeec',
                'epigeec_hdf5' : 'hdf5',
                'epigeec_filtered' : 'filtered',
                'epigeec_output' : 'output',
                'report_dir' : 'report'
                }

        return dirs

    @property
    def contrasts(self):
        contrasts = super(EpiQC, self).contrasts

        # Parse contrasts to retrieve name and type
        for contrast in contrasts:
            if re.search("^\w[\w.-]*,[BN]$", contrast.name):
                contrast.real_name = contrast.name.split(",")[0]
                if contrast.name.split(",")[1] == 'B':
                    contrast.type = 'broad'
                elif contrast.name.split(",")[1] == 'N':
                    contrast.type = 'narrow'
            else:
                _raise(SanitycheckError("Error: contrast name \"" + contrast.name + "\" is invalid (should be <contrast>,B for broad or <contrast>,N for narrow)!"))

        return contrasts

    def createInputInfoFile(self):
        """
            Creates the input info file for ChromImpute after converting the dataset in a bedgraph format.
        """
        with open(config.param('chromimpute', 'inputinfofile'), "w+") as inputinfofile:
            for contrast in self.contrasts:
                if contrast.treatments:
                    for sample in contrast.treatments:
                        inputinfofile.write("{sample}\t{histone}\t{file}".format(
                            sample=sample.name,
                            histone=contrast.real_name,
                            file=os.path.join(self.output_dirs['bedgraph_converted_directory'], sample.name + ".bedgraph")))

        # inputinfofile = open(config.param('chromimpute', 'inputinfofile'), "w+")
        # marks = config.param('chromimpute', 'marks') # Looks for the marks in the ini file (marks are seperated with commas)
        # marks = marks.split(",")
        # cpt = 0

        # # Convert bigwig files to bedgraph and create inputinfofile
        # for sample in self.samples:
        #     for readset in sample.readsets:
        #         if readset.bigwig: # Check if the readset has a BIGWIG column != None
        #             bigwig_file = readset.bigwig
        #             if cpt < len(marks):
        #                 mark = marks[cpt]
        #                 cpt += 1
        #         else:                      # If not, we search for the path from a chipseq pipeline
        #             prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Find path to chipseq folder
        #             bigwig_file = os.path.join(prefix_path, "tracks", sample.name, "bigWig", sample.name + ".bw") # Create path to bigwig file
        #             mark = config.param('DEFAULT', 'chip_type')

        #         # mark = os.path.basename(bigwig_file+".bedgraph").split(".")[-4] # Use only for IHEC database samples
        #         inputinfofile.write(sample.name+"\t"+mark+"\t"+os.path.basename(bigwig_file+".bedgraph")+"\n")

        # inputinfofile.close()

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
            Runs the tool bigWigInfo on bigwig files

            If an EpiQC readset is given to the pipeline, we search for the BIGWIG column to get the files
            If epiqc is ran after a chipseq pipeline, the path to the bigwig files is reconstructed through the location of the chipseq readset file (HAS TO BE IN THE SAME FOLDER AS THE OUTPUT OF THE CHIPSEQ PIPELINE)
        """
        jobs = []

        output_dir = self.output_dirs['bigwiginfo_output_directory']

        log.debug("bigwiginfo_creating_jobs")
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig: # Check if the readset has a BIGWIG column != None
                    bigwig_file = readset.bigwig
                else:                      # If not, we search for the path from a chipseq pipeline
                    prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Find path to chipseq folder
                    bigwig_file = os.path.join(prefix_path, "tracks", sample.name, "bigWig", sample.name + ".bw") # Create path to bigwig file

                file_ext = bigwig_file.split(".")[-1]
                if file_ext not in ("bigWig", "bw"):# != "bigWig" and file_ext != "bw":
                    raise Exception("Error : " + os.path.basename(bigwig_file) + " not a bigWig file !")
                job = concat_jobs([
                    Job(command="mkdir -p " + output_dir),
                    bigwiginfo.bigWigInfo(bigwig_file, self.output_dirs['bigwiginfo_output_directory'])
                    ])
                job.name = "bigwiginfo." + readset.name
                jobs.append(job)

        return jobs

    def bigwig_to_bedgraph(self):
        """
            Converts bigwig files in readset to bedgraph files (Used to train ChromImpute)
        """
        jobs = []

# /!\ gz bedgraph to save space

#         jobs.append(Job(
#             [],
#             ['dataset_dir'],
#             [],
#             command="""\
# mkdir -p \\
#   {dataset_dir} \\
#   {dataset_error}""".format(
#     dataset_dir=self.output_dirs['bedgraph_converted_directory'],
#     dataset_error=self.output_dirs['bedgraph_converted_directory']+"_error"),
#             name='mkdirs_bedgraph_converted_directory'))

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig: # Check if the readset has a BIGWIG column  != None
                    bigwig_file = readset.bigwig
                else:                      # If not, we search for the path from a chipseq pipeline
                    prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Find path to chipseq folder
                    bigwig_file = os.path.join(prefix_path, "tracks", sample.name, "bigWig", sample.name + ".bw") # Create path to bigwig file

            output_bedgraph = os.path.join(self.output_dirs['bedgraph_converted_directory'], sample.name + ".bedgraph")
            output_bedgraph_gz = os.path.join(self.output_dirs['bedgraph_converted_directory'], sample.name + ".bedgraph.gz")

            job = concat_jobs([
                Job(command="mkdir -p " + self.output_dirs['bedgraph_converted_directory']),
                bigwiginfo.bigWigToBedGraph(bigwig_file, output_bedgraph),
                Job(output_files=[output_bedgraph_gz], command="gzip -f " + output_bedgraph)
                ])
            job.name = "bigwig_to_bedgraph." + sample.name
            jobs.append(job)
            # jobs.append(bigwiginfo.bigWigToBedGraph(bigwig_file, self.output_dirs['bedgraph_converted_directory']))

        return jobs

    def chromimpute_preprocess(self):
        """
            Runs the training steps (Convert, ComputeGlobalDist, GenerateTrainData, Train) of the ChromImpute tool on the bigwig files given to the pipeline.
            This step can be skipped if the trained data from the IHEC database is used.
            All the output are stored in the imputation directory.

            Important note :
            The name of the inputinfofile and the dataset, the path to the chromsizes, the resolution and chromosome number has to be specified and the base.ini file under the chromimpute section.

            If the readset file has a BIGWIG column, we use these files for the pipeline
            If epiqc is ran after a chipseq pipeline, the path to the bigwig files is reconstructed through the location of the chipseq readset file (HAS TO BE IN THE SAME FOLDER AS THE OUTPUT OF THE CHIPSEQ PIPELINE)
        """
        jobs = []

        inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'], "inputinfofile.ChromImpute.txt")
        bedgraph_converted_files = []
        for sample in self.samples:
            bedgraph_converted_files.append(os.path.join(self.output_dirs['bedgraph_converted_directory'], sample.name + ".bedgraph.gz"))

        job = concat_jobs([
            Job(input_files=bedgraph_converted_files,
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
            ),
            Job(output_files=[inputinfofile],
                command="""\
cp /dev/null {inputinfofile} && \\
for contrast in {contrasts}
do
  for sample in {samples}
  do
    echo -e "$sample\\t$contrast\\t$sample.bedgraph.gz" >> {inputinfofile}
  done
done""".format(
    inputinfofile=inputinfofile,
    contrasts=" ".join([contrast.real_name for contrast in self.contrasts if contrast.treatments]),
    samples=" ".join([sample.name for contrast in self.contrasts if contrast.treatments for sample in contrast.treatments]),
    path=self.output_dirs['bedgraph_converted_directory'])
            )
        ])

        job.name = "chromimpute_preprocess_train_step"

        jobs.append(job)


        return jobs


    def chromimpute_convert(self):
        """
            Creates a job for chromimpute Convert for each unique mark in the dataset
        """
        jobs = []

        inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'], "inputinfofile.ChromImpute.txt")
        input_dir = self.output_dirs['bedgraph_converted_directory']
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converted_directory'])

        for contrast in self.contrasts:
            if contrast.treatments:
                output_files = []
                for sample in contrast.treatments:
                    with open(os.path.join(os.environ[config.param('chromimpute', 'chromosome_size').split("/")[0].replace("$", "")], config.param('chromimpute', 'chromosome_size').replace(config.param('chromimpute', 'chromosome_size').split("/")[0]+"/", "")), "r") as chrominfofile:
                        for line in chrominfofile:
                            output_files.append(os.path.join(output_dir, "%s_%s.bedgraph.gz.wig.gz" % (line.split("\t")[0], sample.name)))
                        job = chromimpute.convert(input_dir, output_dir, output_files, inputinfofile, contrast.real_name, sample.name)
                        jobs.append(job)
                # treatment_files = [os.path.join(self.output_dirs['bedgraph_converted_directory'], sample.name + ".bedgraph.gz") for sample in contrast.treatments]

        return jobs

    def chromimpute_compute_global_dist(self):
        """
            Creates a job for chromimpute ComputeGlobalDist for each unique mark in the dataset
        """
        jobs = []

        inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'], "inputinfofile.ChromImpute.txt")
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converted_directory'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_distance_directory'])

        for contrast in self.contrasts:
            if contrast.treatments:
                input_files = []
                output_files = []
                for sample in contrast.treatments:
                    with open(os.path.join(os.environ[config.param('chromimpute', 'chromosome_size').split("/")[0].replace("$", "")], config.param('chromimpute', 'chromosome_size').replace(config.param('chromimpute', 'chromosome_size').split("/")[0]+"/", "")), "r") as chrominfofile:
                        for line in chrominfofile:
                            input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converted_directory'], "%s_%s.bedgraph.gz.wig.gz" % (line.split("\t")[0], sample.name)))
                    output_files.append(os.path.join(output_dir, "%s_%s.txt" % (sample.name, contrast.real_name)))
                job = chromimpute.compute_global_dist(input_files, output_dir, output_files, converteddir, inputinfofile, contrast.real_name)
                jobs.append(job)
        # jobs.append(chromimpute.compute_global_dist(input_dir, output_dir, inputinfofile))

        return jobs

    def chromimpute_generate_train_data(self):
        """
            Creates a job for chromimpute GenerateTrainData for each unique mark in the dataset
        """
        jobs = []

        inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'], "inputinfofile.ChromImpute.txt")
        distancedir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_distance_directory'])
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converted_directory'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_traindata_directory'])

        # chroms = []
        # with open(os.path.join(os.environ[config.param('chromimpute', 'chromosome_size').split("/")[0].replace("$", "")], config.param('chromimpute', 'chromosome_size').replace(config.param('chromimpute', 'chromosome_size').split("/")[0]+"/", "")), "r") as chrominfofile:#with open(config.param('chromimpute', 'chrominfofile')) as chrominfofile:             os.environ[config.param('chromimpute', 'chromosome_size').split("/")[0].replace("$", "")]
        #     for line in chrominfofile:
        #         chroms.append(line.split("\t")[0])

        for contrast in self.contrasts:
            if contrast.treatments:
                input_files = []
                output_files = []
                indice = 0
                for sample in contrast.treatments:
                    with open(os.path.join(os.environ[config.param('chromimpute', 'chromosome_size').split("/")[0].replace("$", "")], config.param('chromimpute', 'chromosome_size').replace(config.param('chromimpute', 'chromosome_size').split("/")[0]+"/", "")), "r") as chrominfofile:
                        for line in chrominfofile:
                            input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converted_directory'], "%s_%s.bedgraph.gz.wig.gz" % (line.split("\t")[0], sample.name)))
                    input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_distance_directory'], "%s_%s.txt" % (sample.name, contrast.real_name)))
                    output_files.append(os.path.join(output_dir, "attributes_%s_%i_0.txt.gz" % (contrast.real_name, indice)))
                    output_files.append(os.path.join(output_dir, "traindata_%s_%i_0.txt.gz" % (contrast.real_name, indice)))
                    indice += 1
                # for chrom in chroms:
                jobs.append(chromimpute.generate_train_data(input_files, output_dir, output_files, converteddir, distancedir, inputinfofile, contrast.real_name))

        return jobs

    def chromimpute_train(self):
        """
            Creates a job for chromimpute Train for every sample mark given in the dataset
        """
        jobs = []

        inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'], "inputinfofile.ChromImpute.txt")
        traindatadir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_traindata_directory'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_predictor_directory'])

        for contrast in self.contrasts:
            if contrast.treatments:
                input_files = []
                output_files = []
                indice_i = 0
                for sample in contrast.treatments:
                    output_files = []
                    input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_traindata_directory'], "attributes_%s_%i_0.txt.gz" % (contrast.real_name, indice_i)))
                    input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_traindata_directory'], "traindata_%s_%i_0.txt.gz" % (contrast.real_name, indice_i)))
                    for indice_j, sample_j in enumerate(sorted(contrast.treatments, key=attrgetter("name"))):
                        if sample_j != sample:
                            output_files.append(os.path.join(output_dir, "useattributes_%s_%s_%i_0.txt.gz" % (sample.name, contrast.real_name, indice_j)))
                    indice_i += 1
                    jobs.append(chromimpute.train(input_files, output_dir, output_files, traindatadir, inputinfofile, sample.name, contrast.real_name))

        return jobs

    def chromimpute_train_step(self):
        """
            Runs the training steps (Convert, ComputeGlobalDist, GenerateTrainData, Train) of the ChromImpute tool on the bigwig files given to the pipeline.
            This step can be skipped if the trained data from the IHEC database is used.
            All the output are stored in the imputation directory.

            Important note :
            The name of the inputinfofile and the dataset, the path to the chromsizes, the resolution and chromosome number has to be specified and the base.ini file under the chromimpute section.

            If the readset file has a BIGWIG column, we use these files for the pipeline
            If epiqc is ran after a chipseq pipeline, the path to the bigwig files is reconstructed through the location of the chipseq readset file (HAS TO BE IN THE SAME FOLDER AS THE OUTPUT OF THE CHIPSEQ PIPELINE)
        """
        # TODO idea: Delete predictordir to save space after run

        jobs = []

        # inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'], "inputinfofile.ChromImpute.txt")
        # path_dataset = os.path.join(os.getcwd(), config.param('chromimpute', 'dataset'))

        # self.createInputInfoFile()
        # read_inputinfofile = csv.reader(open(config.param('chromimpute', 'inputinfofile'), 'rb'), delimiter='\t')
        # samplesMarksFiles = self.parseInputInfoFile(read_inputinfofile)

        log.debug("chromimpute_convert")
        jobs.extend(self.chromimpute_convert())
        log.debug("chromimpute_compute_global_dist")
        jobs.extend(self.chromimpute_compute_global_dist())
        log.debug("chromimpute_generate_train_data")
        jobs.extend(self.chromimpute_generate_train_data())
        log.debug("chromimpute_train")
        jobs.extend(self.chromimpute_train())

        return jobs


    def chromimpute_apply(self):
        """
            Creates a job for chromimpute Apply for every sample mark given in the dataset
        """
        jobs = []

        inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'], "inputinfofile.ChromImpute.txt")
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converted_directory'])
        distancedir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_distance_directory'])
        predictordir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_predictor_directory'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_apply'])

        for contrast in self.contrasts:
            if contrast.treatments:
                for sample in contrast.treatments:
                    input_files = []
                    for indice, sample_j in enumerate(sorted(contrast.treatments, key=attrgetter("name"))):
                        if sample_j != sample:
                            input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_predictor_directory'], "useattributes_%s_%s_%i_0.txt.gz" % (sample.name, contrast.real_name, indice)))
                    # output = output_dir+"_"+samplemark[0]+"_"+samplemark[1]
                    jobs.append(chromimpute.apply(input_files, output_dir, converteddir, distancedir, predictordir, inputinfofile, sample.name, contrast.real_name))

        return jobs

    def chromimpute_eval(self, samplesMarksFile):
        """
            Creates a job for chromimpute Eval for every sample mark given in the dataset
        """
        jobs = []

        input_base = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_apply']) # input_base is the name of the folder containing imputed files
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converted_directory'])
        percent1 = config.param('chromimpute', 'percent1')
        percent2 = config.param('chromimpute', 'percent2')

        for contrast in self.contrasts:
            if contrast.treatments:
                input_dir = input_base+"_"+sampleMarkFile[0]+"_"+sampleMarkFile[1]
                output_path = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_eval'], "eval_"+sampleMarkFile[0]+"_"+sampleMarkFile[1]+".txt")
                jobs.append(chromimpute.eval(input_dir, percent1, percent2, converteddir, sampleMarkFile[2]+".wig.gz", input_base, "impute_"+sampleMarkFile[0]+"_"+sampleMarkFile[1]+".wig.gz", output_path))

        return jobs


    def chromimpute_compute_metrics(self):
        """
            Imputes the tracks en runs the Eval step on the bigwig files given to the pipeline.

            All the output are stored in the imputation directory.
            The imputed files can be found in : imputation/output
            The eval files can be found in : imputation/eval

            Important note :
            The name of the inputinfofile and the dataset, the path to the chromsizes, the resolution and chromosome number has to be specified and the base.ini file under the chromimpute section.

            If a readset file is given to the pipeline, we search for the BIGWIG column to get the files
            If epiqc is ran after a chipseq pipeline, the path to the bigwig files is reconstructed through the location of the chipseq readset file (HAS TO BE IN THE SAME FOLDER AS THE OUTPUT OF THE CHIPSEQ PIPELINE)
        """
        jobs = []

#         jobs.append(Job(
#             [],
#             [],
#             [],
#             command="""\
# mkdir -p \\
#   {output_dir}/{apply} \\
#   {output_dir}/{eval}""".format(
#       output_dir=self.output_dirs['chromimpute_output_directory'],
#       apply=self.output_dirs['chromimpute_output'],
#       eval=self.output_dirs['chromimpute_eval']),
#             name='mkdirs_chromimpute_metrics'))

#         path_inputinfofile = os.path.join(os.getcwd(), config.param('chromimpute', 'inputinfofile'))
#         path_dataset = os.path.join(os.getcwd(), config.param('chromimpute', 'dataset'))


        log.debug("chromimpute_apply")
        jobs.extend(self.chromimpute_apply())
        # log.debug("chromimpute_eval")
        # jobs.extend(self.chromimpute_eval())

        return jobs

    def signal_to_noise(self):
        """
            Uses the converted files obtained from chromimpute and calculates the ratio between the top 10% & top 5% of the signals over the sum the the signals
            Outputs the results in a tsv file
        """
        # TODO : Delete decompressed converted files in imputation/converteddir to save space
        jobs = []

        jobs.append(Job(
            [],
            ['signal_noise'],
            [],
            command=
"mkdir \
  {output_dir}".format(
    output_dir=self.output_dirs['signal_to_noise_output_directory']),
            name="mkdir_signal_noise"))

        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converted_directory'])

        log.debug("signal_to_noise")
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig != None: # Check if the readset has a BIGWIG column
                    bigwig_file = readset.bigwig
                    converted_file = os.path.basename(bigwig_file)+".bedgraph.wig.gz"
                else:                      # If not, we search for the path from a chipseq pipeline
                    bigwig_file = sample.name
                    converted_file = os.path.basename(bigwig_file)+".bw.bedgraph.wig.gz"


                output_file = os.path.join(self.output_dirs['signal_to_noise_output_directory'], converted_file+".tsv")

                jobs.append(Job(
                    ['signal_noise', converteddir],
                    [output_file],
                    [],
                    name='signal_noise',
                    command="""\
python ../genpipes/bfx/wigSignalNoise.py \\
  -i {file} \\
  -p1 {percent1} \\
  -p2 {percent2} \\
  -o {output_dir}""".format(
                    file=os.path.join(converteddir, "chr1_" + converted_file),
                    percent1=config.param('signal_noise', 'percent1'),
                    percent2=config.param('signal_noise', 'percent2'),
                    output_dir=output_file
                    )))

        return jobs


    def epigeec_tohdf5(self):
        jobs = []

        input_dir = 'epigeec'
        output_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_hdf5'])

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig != None: # Check if the readset has a BIGWIG column
                    bigwig_file = readset.bigwig
                else:                      # If not, we search for the path from a chipseq pipeline
                    prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Find path to chipseq folder
                    bigwig_file = os.path.join(prefix_path, "tracks", sample.name, "bigWig", sample.name + ".bw") # Create path to bigwig file

                file_ext = bigwig_file.split(".")[-1]
                if file_ext != "bigWig" and file_ext != "bw":
                    raise Exception("Error : " + os.path.basename(bigwig_file) + " not a bigWig file !")

                jobs.append(epigeec.tohdf5(input_dir, output_dir, bigwig_file))

        return jobs

    def epigeec_filter(self, hdf5_files):
        jobs = []

        input_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_hdf5'])
        output_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_filtered'])

        for hdf5_file in hdf5_files:
            path_to_hdf5_file = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_hdf5'], hdf5_file)
            jobs.append(epigeec.filter(input_dir, output_dir, path_to_hdf5_file))

        return jobs

    def epigeec_correlate(self, skip_filter_step, hdf5_files):
        file_list = open(config.param('epigeec', 'hdf5_list'), "w+")

        for hdf5_file in hdf5_files:
            if skip_filter_step:
                path_to_hdf5_file = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_hdf5'], hdf5_file)
            else:
                path_to_hdf5_file = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_filtered'], hdf5_file)

            file_list.write(path_to_hdf5_file + "\n")

        file_list.close()

        if skip_filter_step:
            input_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_hdf5'])
        else:
            input_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_filtered'])
        output_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_output'])

        return epigeec.correlate(input_dir, output_dir, config.param('epigeec', 'hdf5_list'))

    def epigeec(self):
        """
            Runs the epigeec tool on the bigwig files given to the pipeline
            Bigwig files are first converted to the hdf5 format, then the correlation matrix is computed

            If a readset file is given to the pipeline, we search for the BIGWIG column to get the files
            If epiqc is ran after a chipseq pipeline, the path to the bigwig files is reconstructed through the location of the chipseq readset file (HAS TO BE IN THE SAME FOLDER AS THE OUTPUT OF THE CHIPSEQ PIPELINE)
        """
        jobs = []

        skip_filter_step = config.param('epigeec', 'select') == '' and config.param('epigeec', 'exclude') == ''

        filtered_dir = self.output_dirs['epigeec_filtered']
        if skip_filter_step:
            filtered_dir = ""

        jobs.append(Job(
            [],
            ['epigeec'],
            [],
            command=
"mkdir -p \
  {output_dir}/{hdf5_dir} \
  {output_dir}/{filtered_dir} \
  {output_dir}/{correlate}".format(
    output_dir=self.output_dirs['epigeec_output_directory'],
    hdf5_dir=self.output_dirs['epigeec_hdf5'],
    filtered_dir=filtered_dir,
    correlate=self.output_dirs['epigeec_output']),
            name="mkdir_epigeec"))

        hdf5_files = []

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig != None: # Check if the readset has a BIGWIG column
                    bigwig_file = readset.bigwig
                else:                      # If not, we search for the path from a chipseq pipeline
                    prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Find path to chipseq folder
                    bigwig_file = os.path.join(prefix_path, "tracks", sample.name, "bigWig", sample.name + ".bw") # Create path to bigwig file

                hdf5_files.append(os.path.basename(bigwig_file) + ".hdf5")

        log.debug("epigeec_tohdf5")
        jobs.extend(self.epigeec_tohdf5())

        if not skip_filter_step: # We skip this step if there are no filter files
            log.debug("epigee_filter")
            jobs.extend(self.epigeec_filter(hdf5_files))

        log.debug("epigeec_correlate")
        jobs.append(self.epigeec_correlate(skip_filter_step, hdf5_files))

        return jobs

    def epiqc_report(self):
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

        jobs.append(Job(
            [],
            ['epiqc_report'],
            [],
            command=
"mkdir \
  {output_dir}".format(
            output_dir=self.output_dirs['report_dir']),
            name="mkdir_report_epiqc"))

        chromimpute_output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_eval'])

        read_inputinfofile = csv.reader(open(config.param('chromimpute', 'inputinfofile'), 'rb'), delimiter='\t')
        samplesMarksFiles = self.parseInputInfoFile(read_inputinfofile)

        marks = config.param('chromimpute', 'marks')
        marks = marks.split(",")
        cpt = 0

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig is not None:
                    bigwiginfo_file = os.path.join(self.output_dirs['bigwiginfo_output_directory'], "bigwiginfo_"+os.path.basename(readset.bigwig)+".txt")
                    signal_noise_file = os.path.join(self.output_dirs['signal_to_noise_output_directory'], os.path.basename(readset.bigwig)+".bedgraph.wig.gz.tsv")
                    mark = marks[cpt] # Corresponds to the marks specified in the epiqc.base.ini
                    cpt += 1
                else:
                    bigwiginfo_file = os.path.join(self.output_dirs['bigwiginfo_output_directory'], "bigwiginfo_"+sample.name+".bw.txt")
                    signal_noise_file = os.path.join(self.output_dirs['signal_to_noise_output_directory'], sample.name+".bw.bedgraph.wig.gz.tsv")
                    mark = config.param('DEFAULT', 'chip_type')

                report_file = os.path.join(self.output_dirs['report_dir'], "report_"+sample.name+"_"+mark+".txt")

                jobs.append(Job(
                    ['epiqc_report', bigwiginfo_file],
                    [],
                    [],
                    name="report_bigwiginfo_" + os.path.basename(bigwiginfo_file),
                    command=
"python ../genpipes/bfx/epiqc_report.py \
  -b {bigwiginfo_file} \
  -cb {chromCount} \
  -bc1 {low_alert_bases_covered} \
  -bc2 {medium_alert_bases_covered} \
  -o {output_file}".format(
                    bigwiginfo_file=bigwiginfo_file,
                    chromCount=config.param('epiqc_report', 'chromcount_threshold'),
                    low_alert_bases_covered=config.param('epiqc_report', 'low_alert_bases_covered'),
                    medium_alert_bases_covered=config.param('epiqc_report', 'medium_alert_bases_covered'),
                    output_file=report_file)))

                eval_file = "eval_"+sample.name+"_"+mark+".txt"
                eval_file = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_eval'], eval_file)

#                 jobs.append(Job(
#                     ['epiqc_report', eval_file],
#                     [],
#                     [],
#                     name = "report_eval_"+os.path.basename(eval_file),
#                     command =
# "python ../genpipes/bfx/epiqc_report.py \
#   -c {eval_file} \
#   -p1 {percent1} \
#   -p2 {percent2} \
#   -ct1 {chromimpute_threshold_M} \
#   -ct2 {chromimpute_threshold_L} \
#   -o {output_file}".format(
#                     eval_file = eval_file,
#                     percent1 = config.param('chromimpute', 'percent1'),
#                     percent2 = config.param('chromimpute', 'percent2'),
#                     chromimpute_threshold_M=config.param('epiqc_report', 'chromimpute_threshold_M'),
#                     chromimpute_threshold_L=config.param('epiqc_report', 'chromimpute_threshold_L'),
#                     output_file = report_file)))

#             jobs.append(Job(
#                     ['epiqc_report', signal_noise_file],
#                     [],
#                     [],
#                     name = "report_signal_noise_" + os.path.basename(signal_noise_file),
#                     command =
# "python ../genpipes/bfx/epiqc_report.py \
#   -s {signal_noise_file} \
#   -s1 {percent1} \
#   -s2 {percent2} \
#   -st1 {signal_noise_threshold_M} \
#   -st2 {signal_noise_threshold_L} \
#   -o {output_file}".format(
#                     signal_noise_file = signal_noise_file,
#                     percent1 = config.param('signal_noise', 'percent1'),
#                     percent2 = config.param('signal_noise', 'percent2'),
#                     signal_noise_threshold_M = config.param('epiqc_report', 'signal_noise_threshold_M'),
#                     signal_noise_threshold_L = config.param('epiqc_report', 'signal_noise_threshold_L'),
#                     output_file = report_file)))

        input_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_output'])
        jobs.append(Job(
                ['epiqc_report', input_dir],
                [],
                [],
                name="report_heatmap",
                command=
"python ../genpipes/bfx/epiqc_report.py \
  -e {correlation_matrix} \
  -o {output_dir}".format(
                correlation_matrix=input_dir+"/correlation_matrix.tsv",
                output_dir=self.output_dirs['report_dir'])))

        return jobs

    @property
    def steps(self):
        # TODO : - Create steps table
        return [
            self.bigwiginfo,
            self.bigwig_to_bedgraph,
            self.chromimpute_preprocess,
            self.chromimpute_train_step,
            self.chromimpute_compute_metrics,
            # self.signal_to_noise,
            # self.epigeec,
            # self.epiqc_report
            ]

if __name__ == "__main__":
    EpiQC()
