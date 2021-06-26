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
import matplotlib
matplotlib.use('Agg')

import argparse
import csv
import os
import numpy as np
import matplotlib.pyplot as plt



def report_bigwiginfo(bigwiginfo_file, report_file_name):
    report_file = open(report_file_name, "a")

    with open(bigwiginfo_file) as file:
        for line in file:
            line = line.split(":")

            if line[0] == "chromCount":
                chrom_count = int(line[1])
                if chrom_count < 23:
                    report_file.write("BigWigInfo : HIGH LEVEL ALERT : Chromosome count < 23  ! (chromCount = " + str(chrom_count) + "\n")
                    report_file.close()
                    return
                else:
                    report_file.write("BigWigInfo : Chromosome count passed !\n")
                    report_file.close()
                    return

    report_file.write("BigWigInfo : High Level Alert : chromCount not found ! (check if file is .bigwig)\n")
    report_file.close()
    return


def report_chromimpute(chromimpute_eval_file, percent1, percent2, report_file_name):
    report_file = open(report_file_name, "a")

    read_eval_file = csv.DictReader(open(chromimpute_eval_file), delimiter = "\t")

    for line in read_eval_file:
        observed_impute = line['OBSERVED_' + percent1 + '_IMPUTE_' + percent2]
        both = line['BOTH_' + percent1]

        if float(observed_impute) < 30:
            report_file.write("ChromImpute : MEDIUM Level Alert : OBSERVED_"+percent+"_IMPUTE_"+percent2 +"< 30% (OBSERVED_"+percent+"_IMPUTE_"+percent2 + " = " + observed_impute + "%)\n")
        if float(both) < 20:
            report_file.write("ChromImpute : Low Level Alert : BOTH_"+percent1 + "< 20% (BOTH_"+percent1 + " = "+ both + "%)\n")
        else:
            report_file.write("ChromImpute : ChromImpute evaluation passed !\n")

    report_file.close()

def report_signal_noise(signal_noise_file, report_file_name):
    report_file = open(report_file_name, "a")

    read_signal_noise_file = csv.DictReader(open(signal_noise_file), delimiter = "\t")

    for line in read_signal_noise_file:
        ratio_10 = float(line['Ratio top 10% bins'])
        ratio_5 = float(line['Ratio top 5% bins'])

    if ratio_10 < 0.3:
        report_file.write("Signal to noise : MEDIUM Level Alert : signal in top 10% bins < 30% (" + str(ratio_10) +")\n")
    if ratio_5 < 0.2:
        report_file.write("Signal to noise : Low Level Alert : signal in top 5% bins < 20% (" + str(ratio_5) +")\n")
    else:
        report_file.write("Signal to noise : Signal to noise ratio passed ! (" + str(ratio_10) + " > 0.3 & " + str(ratio_5) + " > 0.2)\n")

    report_file.close()

def report_epigeec(correlation_matrix, output_dir):
    file = csv.reader(open(correlation_matrix), delimiter="\t")

    matrix = []
    samples = []
    cpt = 0

    for line in file:
        if cpt == 0:
            samples = line[1:]
            cpt += 1
        else:
            matrix.append(line[1:])

    #Convert matrix to float
    matrix = np.array(matrix)
    matrix = matrix.astype(np.float)

    fig, ax = plt.subplots()
    im = ax.imshow(matrix)

    ax.set_xticks(np.arange(len(samples)))
    ax.set_yticks(np.arange(len(samples)))

    ax.set_xticklabels(samples)
    ax.set_yticklabels(samples)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(samples)):
        for j in range(len(samples)):
            text = ax.text(j, i, matrix[i, j],
                           ha="center", va="center", color="w")

    ax.set_title("Samples correlation matrix")
    fig.tight_layout()
    plt.savefig(output_dir+"/correlation_matrix.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Analyses the output of EpiQC")
    parser.add_argument("-b", "--bigwiginfo", help="Analyses a bigwiginfo file", type=str, default="")
    parser.add_argument("-c", "--chromimpute", help="Analyses a chromimpute eval file", type=str, default="")
    parser.add_argument("-p1", "--percent1", help="Percent1 of chromimpute eval (required if -c)", type=str, default="")
    parser.add_argument("-p2", "--percent2", help="Percent2 of chromimpute eval (required if -c)", type=str, default="")
    parser.add_argument("-s", "--signalnoise", help="Analyses a signal_noise file", type=str, default="")
    parser.add_argument("-e", "--epigeec", help="Creates heatmap from correlation matrix", type=str, default="")
    parser.add_argument("-o", "--output", help="File to write results", type=str, default="")
    args = parser.parse_args()

    report_file = args.output

    if args.bigwiginfo != "":
        if args.output != "":
            report_bigwiginfo(args.bigwiginfo, report_file)
        else:
            print("Specify output file with -o !")
    if args.chromimpute != "":
        if args.percent1 != "" and args.percent2 != "":
            if args.output != "":
                report_chromimpute(args.chromimpute, args.percent1, args.percent2, report_file)
            else:
                print("Specify output file with -o !")
        else:
            print("Missing percent1 and percent2 ! (use -p1 and -p2)")
    if args.signalnoise != "":
        if args.output != "":
            report_signal_noise(args.signalnoise, report_file)
        else:
            print("Specify output file with -o !")
    if args.epigeec != "":
        if args.output != "":
            report_epigeec(args.epigeec, report_file)
        else:
            print("Specify output file with -o !")



