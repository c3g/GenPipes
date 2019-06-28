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

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *

from bfx import bigwiginfo
from bfx import chromimpute
from bfx import wigSignalNoise
from bfx import epigeec
from bfx import epiqc_report
from bfx.readset import *

from pipelines.chipseq import chipseq

log = logging.getLogger(__name__)

class EpiQC(chipseq.ChipSeq):
    """
    TODO: write description of pipeline
    """
    def __init__(self, protocol=None):
        self._protocol = protocol
        super(EpiQC, self).__init__(protocol)

    @property
    def output_dirs(self):
        dirs = {'bigwiginfo_output_directory': 'bigwiginfo',
                'chromimpute_output_directory': 'imputation',
                'bedgraph_dataset' : 'dataset',
                'chromimpute_converteddir' : 'converteddir',
                'chromimpute_distancedir' : 'distancedir',
                'chromimpute_traindatadir' : 'traindatadir',
                'chromimpute_predictordir' : 'predictordir',
                'chromimpute_output' : 'output',
                'chromimpute_eval' : 'eval',
                'signal_to_noise_output_directory': 'signal_to_noise',
                'epigeec_output_directory': 'epigeec',
                'epigeec_hdf5' : 'hdf5',
                'epigeec_filtered' : 'filtered',
                'epigeec_output' : 'output',
                'report_dir' : 'report'
                }

        return dirs

    def parseInputInfoFile(self, inputinfofile):
        """
            Parses a chromimpute input info file
            Each element of the array samplesMarksFiles contains respectivly the name of the sample, its mark, and the name of the file
        """
        samplesMarksFiles = []
        for line in inputinfofile:  
            sampleMarkFile = [line[0], line[1], line[2]]
            samplesMarksFiles.append(sampleMarkFile)

        return samplesMarksFiles

    def bigwiginfo(self):
        """
            Runs the tool bigWigInfo on bigwig files
            Outputs will be stored in a new directory : bigwiginfo

            If an EpiQC readset is given to the pipeline, we search for the BIGWIG column to get the files
            If epiqc is ran after a chipseq pipeline, the path to the bigwig files is reconstructed through the location of the chipseq readset file (HAS TO BE IN THE SAME FOLDER AS THE OUTPUT 
            OF THE CHIPSEQ PIPELINE)
        """
        jobs = []

        jobs.append(Job(
            [],
            ['bigwiginfo_dir'],
            [],
            command = 
"mkdir \
  {output_dir}".format(
    output_dir=self.output_dirs['bigwiginfo_output_directory']),
            name="mkdir_bigwiginfo"))

        log.debug("bigwiginfo_creating_jobs")
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig != None: # Check if the readset has a BIGWIG column 
                    bigwig_path = readset.bigwig
                else:                      # If not, we search for the path from a chipseq pipeline
                    prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Find path to chipseq folder
                    bigwig_path = os.path.join(prefix_path, "tracks", sample.name, "bigWig", sample.name + ".bw") # Create path to bigwig file

                file_ext = bigwig_path.split(".")[-1]
                if file_ext != "bigWig" and file_ext != "bw":
                    raise Exception("Error : " + os.path.basename(bigwig_path) + " not a bigWig file !")

                jobs.append(bigwiginfo.bigWigInfo('bigwiginfo_dir', bigwig_path, self.output_dirs['bigwiginfo_output_directory']))

        return jobs

    def bigwig_to_bedgraph(self):
        """
            Converts bigwig files in readset to bedgraph files (Used to train ChromImpute)
        """
        jobs = []

        jobs.append(Job(
            [],
            ['dataset_dir'],
            [],
            command ="""\
mkdir -p \\
  {dataset_dir}""".format(
    dataset_dir = self.output_dirs['bedgraph_dataset']),
            name = 'mkdirs_bedgraph_dataset'))

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig != None: # Check if the readset has a BIGWIG column 
                    bigwig_path = readset.bigwig
                else:                      # If not, we search for the path from a chipseq pipeline
                    prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Find path to chipseq folder
                    bigwig_path = os.path.join(prefix_path, "tracks", sample.name, "bigWig", sample.name + ".bw") # Create path to bigwig file

                jobs.append(bigwiginfo.bigWigToBedGraph('dataset_dir', bigwig_path, self.output_dirs['bedgraph_dataset']))

        return jobs

    def chromimpute_convert(self, inputinfofile, marks):
        jobs = []

        input_dir = self.output_dirs['bedgraph_dataset']
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converteddir'])

        for mark in marks:
            jobs.append(chromimpute.convert(input_dir, output_dir, inputinfofile, mark))
        
        return jobs

    def chromimpute_compute_global_dist(self, inputinfofile, marks):
        jobs = []

        input_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converteddir'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_distancedir'])

        for mark in marks:
            jobs.append(chromimpute.compute_global_dist(input_dir, output_dir, inputinfofile, mark))

        return jobs

    def chromimpute_generate_train_data(self, inputinfofile, marks):
        jobs = []

        input_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_distancedir'])
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converteddir'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_traindatadir'])
        
        for mark in marks:
            jobs.append(chromimpute.generate_train_data(input_dir, output_dir, converteddir, inputinfofile, mark))

        return jobs

    def chromimpute_train(self, inputinfofile, samplesMarks):
        jobs = []

        input_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_traindatadir'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_predictordir'])
        
        for samplemark in samplesMarks:
            jobs.append(chromimpute.train(input_dir, output_dir, inputinfofile, samplemark[0], samplemark[1]))

        return jobs

    def chromimpute_apply(self, inputinfofile, samplesMarks):
        jobs = []

        input_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_predictordir'])
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converteddir'])
        distancedir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_distancedir'])
        predictordir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_predictordir'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_output'])
        
        for samplemark in samplesMarks:
            jobs.append(chromimpute.apply(input_dir, output_dir, converteddir, distancedir ,predictordir, inputinfofile, samplemark[0], samplemark[1]))
        
        return jobs

    def chromimpute_eval(self, samplesMarksFile):
        jobs = []

        input_dir = output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_output'])
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converteddir'])
        percent1 = config.param('chromimpute', 'percent1')
        percent2 = config.param('chromimpute', 'percent2')

        for sampleMarkFile in samplesMarksFile:
            output_path = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_eval'], "eval_"+sampleMarkFile[0]+"_"+sampleMarkFile[1]+".txt")
            jobs.append(chromimpute.eval(input_dir, percent1, percent2, converteddir, sampleMarkFile[2]+".wig.gz", "impute_"+sampleMarkFile[0]+"_"+sampleMarkFile[1]+".wig.gz", output_path))

        return jobs

    def chromimpute_train_step(self):
        jobs= []

        jobs.append(Job(
            [],
            ['chromimpute_train_dir'],
            [],
            command ="""\
mkdir -p \\
  {output_dir}/{converteddir} \\
  {output_dir}/{compute_global_dist} \\
  {output_dir}/{generate_train_data} \\
  {output_dir}/{train}""".format(
    output_dir = self.output_dirs['chromimpute_output_directory'],
    converteddir = self.output_dirs['chromimpute_converteddir'],
    compute_global_dist = self.output_dirs['chromimpute_distancedir'],
    generate_train_data = self.output_dirs['chromimpute_traindatadir'],
    train = self.output_dirs['chromimpute_predictordir']),
            name = 'mkdirs_chromimpute_train'))

        inputinfofile = open(config.param('chromimpute','inputinfofile'), "w+")
        path_inputinfofile = os.path.join(os.getcwd(), config.param('chromimpute','inputinfofile'))
        path_dataset = os.path.join(os.getcwd(), config.param('chromimpute', 'dataset'))

        marks = config.param('chromimpute', 'marks')
        marks = marks.split(",")
        cpt = 0

        # Convert bigwig files to bedgraph and create inputinfofile
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig != None: # Check if the readset has a BIGWIG column 
                    bigwig_path = readset.bigwig
                    mark = marks[cpt]
                    cpt += 1
                else:                      # If not, we search for the path from a chipseq pipeline
                    prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Find path to chipseq folder
                    bigwig_path = os.path.join(prefix_path, "tracks", sample.name, "bigWig", sample.name + ".bw") # Create path to bigwig file
                    mark = config.param('DEFAULT', 'chip_type')

                inputinfofile.write(sample.name+"\t"+mark+"\t"+os.path.basename(bigwig_path+".bedgraph")+"\n")

        inputinfofile.close()

        read_inputinfofile = csv.reader(open(config.param('chromimpute', 'inputinfofile'), 'rb'), delimiter = '\t')
        samplesMarksFiles = self.parseInputInfoFile(read_inputinfofile)

        unique_marks = []
        for mark in samplesMarksFiles:
            if mark[1] not in unique_marks:
                unique_marks.append(mark[1])
        log.debug("Unique marks : "+ str(unique_marks))

        log.debug("chromimpute_convert")
        jobs.extend(self.chromimpute_convert(path_inputinfofile, unique_marks))
        log.debug("chromimpute_compute_global_dist")
        jobs.extend(self.chromimpute_compute_global_dist(path_inputinfofile, unique_marks))
        log.debug("chromimpute_generate_train_data")
        jobs.extend(self.chromimpute_generate_train_data(path_inputinfofile, unique_marks))
        log.debug("chromimpute_train")
        jobs.extend(self.chromimpute_train(path_inputinfofile, samplesMarksFiles))
        log.debug("chromimpute_apply")

        return jobs

    def chromimpute_compute_metrics(self):
        """
            Runs ChromImpute on the bigwig files given to the pipeline.
            It first converts bigwig files to bedgraph, then executes all the steps of the chromimpute tool (convert to eval)

            All the output are stored in the imputation directory.
            The imputed files can be found in : imputation/output
            The eval files can be found in : imputation/eval

            Important note :
            The name of the inputinfofile and the dataset, the path to the chromsizes, the resolution and chromosome number has to be specified and the base.ini file under the 
            chromimpute section.

            If a readset file is given to the pipeline, we search for the BIGWIG column to get the files
            If epiqc is ran after a chipseq pipeline, the path to the bigwig files is reconstructed through the location of the chipseq readset file (HAS TO BE IN THE SAME FOLDER AS THE OUTPUT 
            OF THE CHIPSEQ PIPELINE)
        """
        # TODO idea: Delete predictordir to save space after run
        jobs= []

        jobs.append(Job(
            [],
            ['chromimpute_metrics_dir'],
            [],
            command ="""\
mkdir -p \\
  {output_dir}/{apply} \\
  {output_dir}/{eval}""".format(
    output_dir = self.output_dirs['chromimpute_output_directory'],
    apply = self.output_dirs['chromimpute_output'],
    eval = self.output_dirs['chromimpute_eval']),
            name = 'mkdirs_chromimpute_metrics'))

        inputinfofile = open(config.param('chromimpute','inputinfofile'), "w+")
        path_inputinfofile = os.path.join(os.getcwd(), config.param('chromimpute','inputinfofile'))
        path_dataset = os.path.join(os.getcwd(), config.param('chromimpute', 'dataset'))

        marks = config.param('chromimpute', 'marks')
        marks = marks.split(",")
        cpt = 0

        # Convert bigwig files to bedgraph and create inputinfofile
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig != None: # Check if the readset has a BIGWIG column 
                    bigwig_path = readset.bigwig
                    mark = marks[cpt]
                    cpt += 1
                else:                      # If not, we search for the path from a chipseq pipeline
                    prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Find path to chipseq folder
                    bigwig_path = os.path.join(prefix_path, "tracks", sample.name, "bigWig", sample.name + ".bw") # Create path to bigwig file
                    mark = config.param('DEFAULT', 'chip_type')

                inputinfofile.write(sample.name+"\t"+mark+"\t"+os.path.basename(bigwig_path+".bedgraph")+"\n")

        inputinfofile.close()

        read_inputinfofile = csv.reader(open(config.param('chromimpute', 'inputinfofile'), 'rb'), delimiter = '\t')
        samplesMarksFiles = self.parseInputInfoFile(read_inputinfofile)

        log.debug("chromimpute_apply")
        jobs.extend(self.chromimpute_apply(path_inputinfofile, samplesMarksFiles))
        log.debug("chromimpute_eval")
        jobs.extend(self.chromimpute_eval(samplesMarksFiles))

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
            command = 
"mkdir \
  {output_dir}".format(
    output_dir = self.output_dirs['signal_to_noise_output_directory']),
            name = "mkdir_signal_noise"))
        
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converteddir'])

        log.debug("signal_to_noise")
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig != None: # Check if the readset has a BIGWIG column 
                    bigwig_path = readset.bigwig
                    converted_file = os.path.basename(bigwig_path)+".bedgraph.wig.gz"
                else:                      # If not, we search for the path from a chipseq pipeline
                    bigwig_path = sample.name
                    converted_file = os.path.basename(bigwig_path)+".bw.bedgraph.wig.gz"

                
                output_file = os.path.join(self.output_dirs['signal_to_noise_output_directory'], converted_file+".tsv")

                jobs.append(Job(
                    ['signal_noise',converteddir],
                    [output_file],
                    [],
                    name = 'signal_noise',
                    command = """\
python ../genpipes/bfx/wigSignalNoise.py \\
  -i {file} \\
  -o {output_dir}""".format(
                    file = os.path.join(converteddir, "chr1_" + converted_file),
                    output_dir = output_file
                    )))

        return jobs


    def epigeec_tohdf5(self):
        jobs = []

        input_dir = 'epigeec'
        output_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_hdf5'])

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig != None: # Check if the readset has a BIGWIG column 
                    bigwig_path = readset.bigwig
                else:                      # If not, we search for the path from a chipseq pipeline
                    prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Find path to chipseq folder
                    bigwig_path = os.path.join(prefix_path, "tracks", sample.name, "bigWig", sample.name + ".bw") # Create path to bigwig file

                file_ext = bigwig_path.split(".")[-1]
                if file_ext != "bigWig" and file_ext != "bw":
                    raise Exception("Error : " + os.path.basename(bigwig_path) + " not a bigWig file !")

                jobs.append(epigeec.tohdf5(input_dir, output_dir, bigwig_path))

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
            If epiqc is ran after a chipseq pipeline, the path to the bigwig files is reconstructed through the location of the chipseq readset file (HAS TO BE IN THE SAME FOLDER AS THE OUTPUT 
            OF THE CHIPSEQ PIPELINE)
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
            command = 
"mkdir -p \
  {output_dir}/{hdf5_dir} \
  {output_dir}/{filtered_dir} \
  {output_dir}/{correlate}".format(
    output_dir = self.output_dirs['epigeec_output_directory'],
    hdf5_dir = self.output_dirs['epigeec_hdf5'],
    filtered_dir = filtered_dir,
    correlate = self.output_dirs['epigeec_output']),
            name="mkdir_epigeec"))

        hdf5_files = []

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig != None: # Check if the readset has a BIGWIG column 
                    bigwig_path = readset.bigwig
                else:                      # If not, we search for the path from a chipseq pipeline
                    prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Find path to chipseq folder
                    bigwig_path = os.path.join(prefix_path, "tracks", sample.name, "bigWig", sample.name + ".bw") # Create path to bigwig file            

                hdf5_files.append(os.path.basename(bigwig_path) + ".hdf5")

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
            command = 
"mkdir \
  {output_dir}".format(
            output_dir = self.output_dirs['report_dir']),
            name = "mkdir_report_epiqc"))

        chromimpute_output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_eval'])

        read_inputinfofile = csv.reader(open(config.param('chromimpute', 'inputinfofile'), 'rb'), delimiter = '\t')
        samplesMarksFiles = self.parseInputInfoFile(read_inputinfofile)

        marks = config.param('chromimpute', 'marks')
        marks = marks.split(",")
        cpt = 0

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig != None:
                    bigwiginfo_file = os.path.join(self.output_dirs['bigwiginfo_output_directory'], "bigwiginfo_"+os.path.basename(readset.bigwig)+".txt")
                    signal_noise_file = os.path.join(self.output_dirs['signal_to_noise_output_directory'], os.path.basename(readset.bigwig)+".bedgraph.wig.gz.tsv")
                    mark = marks[cpt]
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
                    name = "report_bigwiginfo_" + os.path.basename(bigwiginfo_file),
                    command = 
"python ../genpipes/bfx/epiqc_report.py \
  -b {bigwiginfo_file} \
  -o {output_file}".format(
                    bigwiginfo_file = bigwiginfo_file,
                    output_file = report_file)))

                eval_file = "eval_" + sample.name+"_"+mark+".txt"
                eval_file = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_eval'], eval_file)

                jobs.append(Job(
                    ['epiqc_report', eval_file],
                    [],
                    [],
                    name = "report_eval_"+os.path.basename(eval_file),
                    command = 
"python ../genpipes/bfx/epiqc_report.py \
  -c {eval_file} \
  -p1 {percent1} \
  -p2 {percent2} \
  -o {output_file}".format(
                    eval_file = eval_file,
                    percent1 = config.param('chromimpute', 'percent1'),
                    percent2 = config.param('chromimpute', 'percent2'),
                    output_file = report_file)))

                jobs.append(Job(
                        ['epiqc_report', signal_noise_file],
                        [],
                        [],
                        name = "report_signal_noise_" + os.path.basename(signal_noise_file),
                        command = 
"python ../genpipes/bfx/epiqc_report.py \
  -s {signal_noise_file} \
  -o {output_file}".format(
                    signal_noise_file = signal_noise_file,
                    output_file = report_file)))



        return jobs



    @property
    def steps(self):
        # TODO : - Create steps table
        return [
            self.bigwiginfo,
            self.bigwig_to_bedgraph,
            self.chromimpute_train_step,
            self.chromimpute_compute_metrics,
            self.signal_to_noise,
            self.epigeec,
            self.epiqc_report]

if __name__ == "__main__":
    EpiQC()