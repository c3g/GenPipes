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
                'epigeec_output' : 'output'
                }

        return dirs

    def parseInputInfoFile(self, inputinfofile):
        samplesMarksFile = []
        for line in inputinfofile:
            sampleMarkFile = [line[0], line[1], line[2]]
            samplesMarksFile.append(sampleMarkFile)

        return samplesMarksFile

    def bigWigInfo(self):
        jobs = []

        jobs.append(Job(
            [],
            ['bigwiginfo'],
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
                    bigwig_path = os.path.join(prefix_path, "tracks", readset.sample.name, "bigWig", readset.sample.name + ".bw") # Create path to bigwig file

                file_ext = bigwig_path.split(".")[-1]
                if file_ext != "bigWig" and file_ext != "bw":
                    raise Exception("Error : " + os.path.basename(bigwig_path) + " not a bigWig file !")

                jobs.append(bigwiginfo.bigWigInfo(bigwig_path, self.output_dirs['bigwiginfo_output_directory']))

        return jobs

    def chromimpute_convert(self, marks):
        jobs = []

        input_dir = self.output_dirs['chromimpute_output_directory']
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converteddir'])

        for mark in marks:
            jobs.append(chromimpute.convert(input_dir, output_dir, mark))
        
        return jobs

    def chromimpute_compute_global_dist(self, marks):
        jobs = []

        input_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converteddir'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_distancedir'])

        for mark in marks:
            jobs.append(chromimpute.compute_global_dist(input_dir, output_dir, mark))

        return jobs

    def chromimpute_generate_train_data(self, marks):
        jobs = []

        input_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_distancedir'])
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converteddir'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_traindatadir'])
        
        for mark in marks:
            jobs.append(chromimpute.generate_train_data(input_dir, output_dir, converteddir, mark))

        return jobs

    def chromimpute_train(self, samplesMarks):
        jobs = []

        input_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_traindatadir'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_predictordir'])
        
        for samplemark in samplesMarks:
            jobs.append(chromimpute.train(input_dir, output_dir, samplemark[0], samplemark[1]))

        return jobs

    def chromimpute_apply(self, samplesMarks):
        jobs = []

        input_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_predictordir'])
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converteddir'])
        distancedir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_distancedir'])
        predictordir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_predictordir'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_output'])
        
        for samplemark in samplesMarks:
            jobs.append(chromimpute.apply(input_dir, output_dir, converteddir, distancedir ,predictordir, samplemark[0], samplemark[1]))
        
        return jobs

    def chromimpute_eval(self, samplesMarksFile):
        jobs = []

        input_dir = output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_output'])
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converteddir'])
        percent1 = config.param('chromimpute', 'percent1')
        percent2 = config.param('chromimpute', 'percent2')

        for sampleMarkFile in samplesMarksFile:
            output_path = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_eval'], "eval_"+sampleMarkFile[0]+"_"+sampleMarkFile[1])
            jobs.append(chromimpute.eval(input_dir, percent1, percent2, converteddir, sampleMarkFile[2]+".wig.gz", "impute_"+sampleMarkFile[0]+"_"+sampleMarkFile[1]+".wig.gz", output_path))

        return jobs

    def chromimpute(self):
        # TODO idea: Automatically convert bigwig to bedgraph, create inputinfofile, delete predictordir to save space
        # Note : path to inputinfofile, dataset and chromsizes, resolution sizes and chrom number are in .ini file
        jobs= []

        jobs.append(Job(
            [],
            [self.output_dirs['chromimpute_output_directory']],
            [],
            command ="""\
mkdir -p \
  {output_dir}/{converteddir} \\
  {output_dir}/{compute_global_dist} \\
  {output_dir}/{generate_train_data} \\
  {output_dir}/{train} \\
  {output_dir}/{apply} \\
  {output_dir}/{eval}""".format(
    output_dir = self.output_dirs['chromimpute_output_directory'],
    converteddir = self.output_dirs['chromimpute_converteddir'],
    compute_global_dist = self.output_dirs['chromimpute_distancedir'],
    generate_train_data = self.output_dirs['chromimpute_traindatadir'],
    train = self.output_dirs['chromimpute_predictordir'],
    apply = self.output_dirs['chromimpute_output'],
    eval = self.output_dirs['chromimpute_eval']),
            name = 'mkdirs_chromimpute'))

        read_inputinfofile = csv.reader(open(config.param('chromimpute', 'inputinfofile'), 'rb'), delimiter = '\t')
        samplesMarksFiles = self.parseInputInfoFile(read_inputinfofile)
        
        unique_marks = []
        for mark in samplesMarksFiles:
            if mark[1] not in unique_marks:
                unique_marks.append(mark[1])
        log.debug("Unique marks : "+ str(unique_marks))

        log.debug("chromimpute_convert")
        jobs.extend(self.chromimpute_convert(unique_marks))
        log.debug("chromimpute_compute_global_dist")
        jobs.extend(self.chromimpute_compute_global_dist(unique_marks))
        log.debug("chromimpute_generate_train_data")
        jobs.extend(self.chromimpute_generate_train_data(unique_marks))
        log.debug("chromimpute_train")
        jobs.extend(self.chromimpute_train(samplesMarksFiles))
        log.debug("chromimpute_apply")
        jobs.extend(self.chromimpute_apply(samplesMarksFiles))
        log.debug("chromimpute_eval")
        jobs.extend(self.chromimpute_eval(samplesMarksFiles))

        return jobs

    def signal_to_noise(self):
        # TODO : Delete decompressed converted files in imputation/converteddir to save space
        # Use converted files from chromimpute to calcule signal to noise ratio
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
        for file in os.listdir(config.param('chromimpute', 'dataset')): # Create path to files that are going to be created by chromimpute_convert
            jobs.append(Job(
                ['signal_noise',converteddir],
                ['signal_noise'],
                [],
                name = 'signal_noise',
                command = """\
python ../genpipes/bfx/wigSignalNoise.py \\
  -i {file} \\
  -o {output_dir}""".format(
                file = os.path.join(converteddir, file + ".wig.gz"),
                output_dir = os.path.join(self.output_dirs['signal_to_noise_output_directory'], file + ".txt")
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
                    bigwig_path = os.path.join(prefix_path, "tracks", readset.sample.name, "bigWig", readset.sample.name + ".bw") # Create path to bigwig file

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
        file_list = open("epigeec_list.txt", "w+")

        for hdf5_file in hdf5_files:
            if skip_filter_step:
                path_to_hdf5_file = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_hdf5'], hdf5_file)
            else:
                path_to_hdf5_file = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_filtered'], hdf5_file)

            file_list.write(path_to_hdf5_file + "\n")

        if skip_filter_step:
            input_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_hdf5'])
        else:
            input_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_filtered'])
        output_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_output'])

        return epigeec.correlate(input_dir, output_dir, "epigeec_list.txt")

    def epigeec(self):
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
                    bigwig_path = os.path.join(prefix_path, "tracks", readset.sample.name, "bigWig", readset.sample.name + ".bw") # Create path to bigwig file            

                hdf5_files.append(os.path.basename(bigwig_path) + ".hdf5")

        jobs.extend(self.epigeec_tohdf5())

        if not skip_filter_step: # We skip this step if there are no filter files
            jobs.extend(self.epigeec_filter(hdf5_files))

        jobs.append(self.epigeec_correlate(skip_filter_step, hdf5_files))

        return jobs

    @property
    def steps(self):
        # TODO : - Create steps table
        return [
            self.bigWigInfo,
            self.chromimpute,
            self.signal_to_noise,
            self.epigeec]

if __name__ == "__main__":
    EpiQC()