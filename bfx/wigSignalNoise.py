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

import gzip
import re
import argparse
import os


def getFileExtension(file):
    parsedFile = file.split(".")
    return parsedFile[-1]


def decompressFile(file):
    decompressedFile = gzip.GzipFile(file, 'rb')
    content = decompressedFile.read()
    decompressedFile.close()

    wigFile = open(file[:-3], "w+")
    wigFile.write(content)
    wigFile.close()
    return filename[:-3]


def sortSignals(filename):
    """
        Sorts the signals of a wig file from strongest to weakest
    """

    # if file gziped we create a new decompressed file
    fileExt = getFileExtension(filename)
    if fileExt == "gz":
        filename = decompressFile(filename)
    fileExt = getFileExtension(filename)
    if not (fileExt == "wig" or fileExt == "bedgraph"):
        raise Exception("Error : " + filename + " has to be .wig or .bedgraph")

    lines = []  # An array containing each lines of the file

    with open(filename) as file:
        for line in file:
            parsed_line = line.split("\t")
            if not re.search("[a-zA-Z]", parsed_line[
                -1]):  # checks if line is a description line
                lines.append(parsed_line[-1])

    lines.sort(reverse=True)
    return lines


def computeMetrics(sorted_signals, percent1, percent2):
    """
        Calculates the sum of the signals and the sum of the top 10 and 5% signals
    """

    signalSum = 0
    sumTopBins_1 = 0
    sumTopBins_2 = 0

    # Creates all the values needed to calculate metrics
    nbBins = len(sorted_signals)
    nbTopBins_1 = int(round(nbBins * percent1))
    nbTopBins_2 = int(round(nbBins * percent2))

    for i in range(nbBins):
        signalSum += float(sorted_signals[i])

    for i in range(nbTopBins_1):
        if (i == nbTopBins_2):
            sumTopBins_2 = sumTopBins_1
        sumTopBins_1 += float(sorted_signals[i])

    return [signalSum, sumTopBins_1, sumTopBins_2]


def outputFile(signalSum, sumTopBins_1, sumTopBins_2, percent1, percent2, output_path):
    """
        Stores the ratio between the signal sum and the top 10 and 5% in an output file
    """
    percent1 = percent1 * 100
    percent2 = percent2 * 100

    columns = "Signal sum\tTop " + str(percent1) + "% bins sum\tTop " + str(percent2) + "% bins sum\tRatio top " + str(
        percent1) + "% bins\tRatio top " + str(percent2) + "% bins"
    if signalSum != 0:
        values = str(signalSum) + "\t" + str(sumTopBins_1) + "\t" + str(sumTopBins_2) + "\t" + str(
            sumTopBins_1 / signalSum) + "\t" + str(sumTopBins_2 / signalSum)
    else:
        values = "0\t0\t0\t0\t0"

    complete_file_name = output_path

    # Writes metrics in a new file
    output = open(complete_file_name, "w+")
    output.write(columns + "\n")
    output.write(values)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description="Calculates the amount of noise in a signal track")
    parser.add_argument("-i", "--input_file", help="Path to input file to analyse", type=str, required=True)
    parser.add_argument("-p1", "--percent1", help="Percent 1 for ratio1", type=float, required=True)
    parser.add_argument("-p2", "--percent2", help="Percent 2 for ratio2", type=float, required=True)
    parser.add_argument("-o", "--output", help="Output directory", type=str, required=False, default=".")
    args = parser.parse_args()

    filename = args.input_file  # converted file
    output_path = args.output
    percent1 = args.percent1
    percent2 = args.percent2

    sorted_signals = sortSignals(filename)
    metrics = computeMetrics(sorted_signals, percent1, percent2)
    outputFile(metrics[0], metrics[1], metrics[2], percent1, percent2, output_path)

