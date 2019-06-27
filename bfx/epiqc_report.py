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

import os
import csv
import argparse

from core.config import *

def report_bigwiginfo(bigwiginfo_file, report_file):
    report_file = open(report_file, "w+")

    with open(bigwiginfo_file) as file:
        for line in file:
            line = line.split(":")

            if line[0] == "chromCount":
                chrom_count = int(line[1])
                if chrom_count < 23:
                    report_file.write("High Level Alert : Chromosome count < 23 (chromCount = " + str(chrom_count) + ")")
                else:
                    report_file.write("Chromosome count passed !")

    report_file.write("chromCount not found !")

def report_chromimpute(chromimpute_eval_file, report_file):
    report_file = open(report_file, "w+")
    read_eval_file = csv.DictReader(open(chromimpute_eval_file), delimiter = "\t")

    percent1 = config.param('chromimpute', 'percent1')
    percent2 = config.param('chromimpute', 'percent2')

    for line in read_eval_file:
        observed_impute = line['OBSERVED_' + percent1 + '_IMPUTE_' + percent2]
        both = line['BOTH_' + percent1]

        if observed_impute < 30:
            report_file.write("Medium Level Alert : OBSERVED_"+percent1+"_IMPUTE_"+percent2 +"< 30% (" + observed_impute + ")")
        if both < 20:
            report_file.write("Low Level Alert : BOTH_"+percent1 + "< 20% (" + both + ")")
        else:
            report_file.write("ChromImpute evaluation passed !")

def report_signal_noise(signal_noise_file, report_file):
    read_signal_noise_file = csv.DictReader(open(signal_noise_file), delimiter = "\t")

    for line in read_signal_noise_file:
        ratio_10 = float(line['Ratio top 10% bins'])
        ratio_5 = float(line['Ratio top 5% bins'])

    if ratio_10 < 0.3:
        report_file.write("Medium Level Alert : signal in top 10% bins < 30% (" + str(ratio_10) +")")
    if ratio_5 < 0.2:
        report_file.write("Low Level Alert : signal in top 5% bins < 20% (" + str(ratio_5) +")")
    else:
        report_file.write("Signal to noise ratio passed !")

def report_epigeec(correlation_matrix):
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Analyses the output of EpiQC")
    parser.add_argument("-b", "--bigwiginfo", help="Analyses a bigwiginfo file", type=str, default="")
    parser.add_argument("-c", "--chromimpute", help="Analyses a chromimpute eval file", type=str, default="")
    parser.add_argument("-s", "--signalnoise", help="Analyses a signal_noise file", type=str, default="")
    parser.add_argument("-e", "--epigeec", help="Creates heatmap from correlation matrix", type=str, default="")
    parser.add_argument("-o", "--output", help="File to write results", type=str, required=True)
    args = parser.parse_args()

    report_file = args.output

    if args.bigwiginfo != "":
        report_bigwiginfo(args.bigwiginfo, report_file)
    if args.chromimpute != "":
        report_chromimpute(args.chromimpute, report_file)
    if args.signalnoise != "":
        report_signal_noise(args.signalnoise, report_file)
    if args.epigeec != "":
        report_epigeec(args.epigeec)




