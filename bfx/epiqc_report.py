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



def report_bigwiginfo(bigwiginfo_file, chromCount, low_alert_bases_covered, medium_alert_bases_covered, report_file_name):
    """
        Reads the output from the BigWigInfo step of the pipeline and determines if there is a highlevel alert or not.
        Chromosome count < 23 || no chromcount found => High level alert
    """
    report_file = open(report_file_name, "a")
    try:
        with open(bigwiginfo_file) as file:
            for line in file:

                line = line.split(":")
                if len(line) == 1:
                    report_file.write("BigWigInfo : HIGH LEVEL ALERT : File might not be a bigwig file !!\n")
                    return

                if line[0] == "chromCount":
                    chrom_count = int(line[1])
                    if chrom_count < int(chromCount):
                        msg = "BigWigInfo : HIGH LEVEL ALERT : Chromosome count < {chromCount}  ! (chromCount = {chromCount})\n".format(
                                chromCount = chromCount
                                )
                        report_file.write(msg)
                    else:
                        report_file.write("BigWigInfo : Chromosome count passed !\n")

                if line[0] == "basesCovered":
                    bases_covered = int(line[1].replace(",",""))
                    if bases_covered < int(low_alert_bases_covered):
                        msg = "BigWigInfo : LOW LEVEL ALERT : Bases covered < {low_alert_bases_covered} ! (bases covered = {bases_covered})\n".format(
                                low_alert_bases_covered = low_alert_bases_covered,
                                bases_covered = bases_covered
                                )
                        report_file.write(msg)
                    elif bases_covered < int(medium_alert_bases_covered):
                        msg = "BigWigInfo : MEDIUM LEVEL ALERT : Bases covered < {medium_alert_bases_covered} ! (bases covered = {bases_covered})\n".format(
                                medium_alert_bases_covered = medium_alert_bases_covered,
                                bases_covered = bases_covered
                                )
                        report_file.write(msg)
                    else:
                        report_file.write("BigWigInfo : Bases covered normal !\n")
    except:
        report_file.write("File " + bigwiginfo_file + " could not be read !")

    report_file.close()
    return


def report_chromimpute(chromimpute_eval_file, percent1, percent2, thresholdM, thresholdL, report_file_name):
    """
        Reads the eval files obtained with ChromImpute and determines if there is a medium/low level alert or not.
        OBSERVED_1_IMPUTE_5 < 30 => Medium level alert
        BOTH_1 < 20 => Low level alert
    """
    report_file = open(report_file_name, "a")

    try:
        read_eval_file = csv.DictReader(open(chromimpute_eval_file), delimiter = "\t")
    except:
        report_file.write("File " + chromimpute_eval_file + " could not be read !")

    for line in read_eval_file:
        observed_impute = float(line['OBSERVED_' + percent1 + '_IMPUTE_' + percent2])
        both = float(line['BOTH_' + percent1])

        if observed_impute < float(thresholdM):
            msg = "ChromImpute : MEDIUM Level Alert : OBSERVED_{percent1}_IMPUTE_{percent2} < {thresholdM}% (OBSERVED_{percent1}_IMPUTE_{percent2} = {observed_impute}%)\n".format(
                    percent1 = percent1,
                    percent2 = percent2,
                    thresholdM = thresholdM,
                    observed_impute = observed_impute
                    )
            report_file.write(msg)
        elif both < float(thresholdL):
            msg = "ChromImpute : Low Level Alert : BOTH_{percent1} < {thresholdL}% (BOTH_{percent1} = {both}%)\n".format(
                    percent1 = percent1,
                    percent2 = percent2,
                    thresholdL = thresholdL,
                    both = both                    
                    )
            report_file.write(msg)
        else:
            msg = "ChromImpute : ChromImpute evaluation passed ! ({observed_impute} > {thresholdM}% & {both} > {thresholdL}%)\n".format(
                    observed_impute = observed_impute,
                    both = both,
                    thresholdM = thresholdM,
                    thresholdL = thresholdL
                    )
            report_file.write(msg)

    report_file.close()

def report_signal_noise(signal_noise_file, toppercent1, toppercent2, thresholdM, thresholdL, report_file_name):
    """
        Reads the files created in the signal to noise step and determines what level alert the file is.
        signal in top 10% bins < 30% => Medium level alert
        signal in top 5% bins < 20% => Low level alert

    """
    report_file = open(report_file_name, "a")

    try:
        read_signal_noise_file = csv.DictReader(open(signal_noise_file), delimiter = "\t")
    except:
        report_file.write("File " + signal_noise_file + " could not be read !")

    toppercent1 = float(toppercent1)*100
    toppercent2 = float(toppercent2)*100
    thresholdL_percent = float(thresholdL)*100
    thresholdM_percent = float(thresholdM)*100

    for line in read_signal_noise_file:
        ratio_1 = float(line['Ratio top '+ str(toppercent1) +'% bins'])
        ratio_2 = float(line['Ratio top '+ str(toppercent2) +'% bins'])

    if ratio_1 < float(thresholdM):
        msg = "Signal to noise : MEDIUM Level Alert : signal in top {toppercent1}% bins < {thresholdM_percent}% ({ratio_1})\n".format(
                toppercent1 = toppercent1,
                thresholdM_percent = thresholdM_percent,
                ratio_1 = ratio_1
                )
        report_file.write(msg)
    elif ratio_2 < float(thresholdL):
        msg = "Signal to noise : Low Level Alert : signal in top {toppercent2}% bins < {thresholdL_percent}% ({ratio_2})\n".format(
                toppercent2 = toppercent2,
                thresholdL_percent = thresholdL_percent,
                ratio_2 = ratio_2
                )
        report_file.write(msg)
    else:
        msg = "Signal to noise : Signal to noise ratio passed ! ( {ratio_1} > {thresholdM_percent} & {ratio_2} > {thresholdL_percent})\n".format(
                ratio_1 = ratio_1,
                thresholdM_percent = thresholdM_percent,
                ratio_2 = ratio_2,
                thresholdL_percent = thresholdL_percent
                )
        report_file.write(msg)

    report_file.close()

def report_epigeec(correlation_matrix, output_dir):
    """
        Creates a heatmap from the correlation matrix obtained with epigeec.
    """
    file = csv.reader(open(correlation_matrix), delimiter="\t")

    matrix = []
    samples = []
    cpt = 0

    for line in file:
        if cpt == 0:
            samples = line[1:]
            cpt = 1
        else:
            matrix.append(line[1:])

    cpt = 0
    for sample in samples:
        split = sample.split(".")
        samples[cpt] = split[2]+"_"+split[-3]
        cpt += 1

    #Convert matrix to float
    matrix = np.array(matrix)
    matrix = matrix.astype(np.float)
    
    cell_size_inch = 0.3
    matrix_height = cell_size_inch * matrix.shape[0]

    fig, ax = plt.subplots()
    if matrix_height > 10:
      fig.set_size_inches(matrix_height, matrix_height+5)
    ax.matshow(matrix)
    ax.set(xticks=np.arange(len(samples)), xticklabels=samples,
           yticks=np.arange(len(samples)), yticklabels=samples)

    plt.setp(ax.get_xticklabels(), rotation=-45, ha="right",
         rotation_mode="anchor")

    fig.tight_layout()
    plt.savefig(output_dir+"/correlation_matrix.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Analyses the output of EpiQC")
    parser.add_argument("-b", "--bigwiginfo", help="Analyses a bigwiginfo file", type=str, default="")
    parser.add_argument("-cb", "--chromcount", help="Chromosome count threshold", type=str, default="23")
    parser.add_argument("-bc1", "--basecovertedL", help="Low level threshold for bases covered (required if -b)", type=str, default="75000000")
    parser.add_argument("-bc2", "--basecovertedM", help="Medium level threshold for bases covered (required if -b)", type=str, default="25000000")
    parser.add_argument("-c", "--chromimpute", help="Analyses a chromimpute eval file", type=str, default="")
    parser.add_argument("-p1", "--percent1", help="Percent1 of chromimpute eval (required if -c)", type=str, default="")
    parser.add_argument("-p2", "--percent2", help="Percent2 of chromimpute eval (required if -c)", type=str, default="")
    parser.add_argument("-ct1", "--cthresholdM", help="Chromimpute threshold medium alert", type=str, default="30")
    parser.add_argument("-ct2", "--cthresholdL", help="Chromimpute threshold low alert", type=str, default="20")
    parser.add_argument("-s", "--signalnoise", help="Analyses a signal_noise file", type=str, default="")
    parser.add_argument("-s1", "--toppercent1", help="Top percent1 for signal to noise", type=str, default="0.1")
    parser.add_argument("-s2", "--toppercent2", help="Top percent2 for signal to noise", type=str, default="0.05")
    parser.add_argument("-st1", "--sthresholdM", help="Signal noise threshold medium alert", type=str, default="0.3")
    parser.add_argument("-st2", "--sthresholdL", help="Signal noise threshold low alert", type=str, default="0.2")
    parser.add_argument("-e", "--epigeec", help="Creates heatmap from the epigeec correlation matrix", type=str, default="")
    parser.add_argument("-o", "--output", help="File to write results", type=str, default="")
    args = parser.parse_args()

    report_file = args.output

    if args.bigwiginfo != "":
        if args.output != "":
            report_bigwiginfo(args.bigwiginfo, args.chromcount, args.basecovertedL, args.basecovertedM, report_file)
        else:
            print("Specify output file with -o !")
    if args.chromimpute != "":
        if args.percent1 != "" and args.percent2 != "":
            if args.output != "":
                report_chromimpute(args.chromimpute, args.percent1, args.percent2, args.cthresholdM, args.cthresholdL, report_file)
            else:
                print("Specify output file with -o !")
        else:
            print("Missing percent1 and percent2 ! (use -p1 and -p2)")
    if args.signalnoise != "":
        if args.toppercent1 != "" and args.toppercent2 != "":
            if args.output != "":
                report_signal_noise(args.signalnoise, args.toppercent1, args.toppercent2, args.sthresholdM, args.sthresholdL, report_file)
            else:
                print("Specify output file with -o !")
        else:
            print("Missing top percents for signal to noise !")
    if args.epigeec != "":
        if args.output != "":
            report_epigeec(args.epigeec, report_file)
        else:
            print("Specify output file with -o !")



