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

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *

from bfx import bigwiginfo
from bfx import chromimpute
from bfx import signal_noise
from bfx import epigeec
from bfx.readset import *

from pipelines import common

log = logging.getLogger(__name__)

class EpiQC(common.Illumina):
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
                'signal_to_noise_output_directory': 'signal_to_noise',
                'epiGeEc_output_directory': 'epiGeEC'
                }

        return dirs

    def bigWigInfo(self):
        jobs = []

        jobs.append(Job(command = "mkdir {output_dir}".format(output_dir=self.output_dirs['bigwiginfo_output_directory']), name="mkdir_bigwiginfo"))

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig != None: # Check if the readset has a BIGWIG column 
                    bigwig_path = readset.bigwig
                else:                       # If not, we search for the path from a chipseq pipeline
                    prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Get path to chipseq folder
                    bigwig_path = os.path.join(prefix_path, "tracks", readset.sample.name, "bigWig", readset.sample.name + ".bw") # Get path to bigwig file

                file_ext = bigwig_path.split(".")[-1]
                if file_ext != "bigWig" and file_ext != "bw":
                    raise Exception("Error : " + os.path.basename(bigwig_path) + " not a bigWig file !")

                jobs.append(bigwiginfo.bigWigInfo(bigwig_path, self.output_dirs['bigwiginfo_output_directory']))

        return jobs

    def chromimpute_convert(self):
        # TODO: - Create inputinfofile.txt with createInputInfo.py
        
        jobs = []
        
        jobs.append(Job(command = "mkdir {output_dir}".format(output_dir=self.output_dirs['chromimpute_output_directory']), name="mkdir_chromimpute"))

        chromimpute.createInputInfo(self.samples, ".")

    def chromimpute_compute_global_dist(self, converteddir, inputInfo):
        jobs = []

    def chromimpute_generate_train_data(self, converteddir, distancedir, inputInfo, mark):
        jobs = []

    def chromimpute_train(self, traindatadir, inputInfo, sample, mark):
        jobs = []

    def chromimpute_apply(self, converteddir, distancedir, predictordir, inputInfo, sample, mark):
        jobs = []

    def chromimpute_eval(self, converteddir, convertedFile, imputedir, imputeFile):
        jobs = []

    def chromimpute(self):
        jobs= []

        self.chromimpute_convert()

    def signal_to_noise(self, signalFile):
        # TODO: - Convert bigWig to wig if not already done
        #       - Create command to run signal_noise.py || Transform signal_noise.py to launch here
        jobs = []

    def epigeec_tohdf5(self, dataset, resolution):
        # TODO: - Create command to convert bigWig to hdf5
        #       - Create txt file containing path to each hdf5 files (one per line)
        #       - Create command to correlate hdf5 files
        jobs = []

    # def epigeec_filter():

    def epigeec_correlate(self, chromSizes):
        jobs = []

    def epigeec_average_score(self, dataset, chromSizes, resolution):
        jobs = []

    @property
    def steps(self):
        # TODO : - Create steps table
        return [
            self.bigWigInfo,
            self.chromimpute,
            self.signal_to_noise,
            self.epigeec_average_score]

if __name__ == "__main__":
    EpiQC()