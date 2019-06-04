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
import os
import sys
import logging

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *

from bfx import bigwiginfo
from bfx import chromimpute
from bfx import signal_noise
from bfx import epigeec

log = logging.getLogger(__name__)

class EpiQC(Pipeline):
    """
    TODO: write description of pipeline
    """
    def __init__(self, protocol=None):
        self._protocol = protocol
        self.version = 1.0
        self.argparser.add_argument("--dataset", help="List of path to signal files (one per line)", type=str, required=True)
        super(EpiQC, self).__init__()

    @property
    def output_dirs(self):
        dirs = {'bigwiginfo_output_directory': 'bigwiginfo',
                'chromimpute_output_directory': 'imputation',
                'signal_to_noise_output_directory': 'signal_to_noise',
                'epiGeEc_output_directory': 'epiGeEC'
                }

        return dirs

    def bigWigInfo(self):
        #TODO: - Convert wig to bigWig if not already done
        #      - Create command for bigWigInfo, output in text file
        dataset = self.args.dataset
        jobs = []

        with open(dataset) as list_signal_file:
            for line in list_signal_file:
                if line.split(".")[-1].rstrip("\n") == "wig":
                    converToBigWig = bigwiginfo.wigToBigWig(line, config.param("epigeec", "chromsizes"))
                    line = line + ".bigWig"
                    jobs.append(converToBigWig)

                if line.split(".")[-1].rstrip("\n") != "bigWig" and line.split(".")[-1].rstrip("\n") != "bw":
                    raise Exception("bigWigInfo : Not a bigWig file !")

                jobs.append(bigwiginfo.bigWigInfo(line))

        return jobs

    def chromimpute_convert(self):
        #TODO: - Create inputinfofile.txt with createInputInfo.py
        dataset = self.args.dataset
        jobs = []

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

    def chromimpute(self, dataset):
        jobs= []

    def signal_to_noise(self, signalFile):
        #TODO: - Convert bigWig to wig if not already done
        #      - Create command to run signal_noise.py || Transform signal_noise.py to launch here
        jobs = []

    def epigeec_tohdf5(self, dataset, resolution):
        #TODO: - Create command to convert bigWig to hdf5
        #      - Create txt file containing path to each hdf5 files (one per line)
        #      - Create command to correlate hdf5 files
        jobs = []

    # def epigeec_filter():

    def epigeec_correlate(self, chromSizes):
        jobs = []

    def epigeec_average_score(self, dataset, chromSizes, resolution):
        jobs = []

    @property
    def steps(self):
        #TODO : - Create steps table
        return [
            self.bigWigInfo,
            self.chromimpute,
            self.signal_to_noise,
            self.epigeec_average_score]

if __name__ == "__main__":
    EpiQC()