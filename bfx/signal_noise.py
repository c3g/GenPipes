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
    # split the file name by .
    parsedFile = file.split(".")
    # return the last item in the list
    return parsedFile[-1]


def decompressFile(file):
    decompressedFile = gzip.GzipFile(file, 'rb')
    content = decompressedFile.read()
    decompressedFile.close()

    #:-3 gives the file name without last 3 digits
    #print(os.path.basename(file))
    wigFile = open(file[:-3], "wb+")
    #print(wigFile)
    # write contents in decompressed wig file into a new file
    wigFile.write(content)
    wigFile.close()
    return filename[:-3]


def removeDecompressFile(filename):
    # get the file extension
    fileExt = getFileExtension(filename)
    # if the file is a gz, is has already been unzipped by this python script. Therefore, we need to remove it.
    # If the file extension is gz then we can remove the uncompressed file
    if fileExt == "gz":
        if os.path.exists(filename[:-3]):
            os.remove(filename[:-3])


def sortSignals(filename):
    # if file gziped we create a new decompressed file
    fileExt = getFileExtension(filename)
    if fileExt == "gz":
        filename = decompressFile(filename)
    fileExt = getFileExtension(filename)
    if not (fileExt == "wig" or fileExt == "bedgraph" or fileExt == "bdg"):
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


def computeMetrics(sorted_signals, perc1, perc2):

    perc1 = float(float(perc1)/100)
    perc2 = float(float(perc2) / 100)
    signalSum = 0
    sumTopBins_10 = 0
    sumTopBins_05 = 0

    # Creates all the values needed to calculate metrics
    nbBins = len(sorted_signals)
    nbTopBins_10 = int(round(nbBins * (perc1)))
    nbTopBins_05 = int(round(nbBins * (perc2)))
    for i in range(nbBins):
        signalSum += float(sorted_signals[i])

    for i in range(nbTopBins_10):
        if (i == nbTopBins_05):
            sumTopBins_05 = sumTopBins_10
        sumTopBins_10 += float(sorted_signals[i])
    return [signalSum, sumTopBins_10, sumTopBins_05]


def outputFile(signalSum, sumTopBins_10, sumTopBins_05, output_path):
    columns = "Signal sum\tTop 10% bins sum\tTop 5% bins sum\tRatio top 10% bins\tRatio top 5% bins"
    values = str(signalSum) + "\t" + str(sumTopBins_10) + "\t" + str(sumTopBins_05) + "\t" + str(
        sumTopBins_10 / signalSum) + "\t" + str(sumTopBins_05 / signalSum)

    complete_file_name = output_path

    # Writes metrics in a new file
    output = open(complete_file_name, "w+")
    output.write(columns + "\n")
    output.write(values)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description="Calculates the amount of noise in a signal track")
    parser.add_argument("-i", "--input_file", help="Path to input file to analyse", type=str, required=True)
    parser.add_argument("-o", "--output", help="Output directory", type=str, required=False, default=".")
    parser.add_argument("-p1", "--percent1", help="Percentage 1", type=str, required=False, default=".")
    parser.add_argument("-p2", "--percent2", help="Percentage 2", type=str, required=False, default=".")
    args = parser.parse_args()

    filename = args.input_file  # converted file
    output_path = args.output
    percent1 = args.percent1
    percent2 = args.percent2
    sorted_signals = sortSignals(filename)
    metrics = computeMetrics(sorted_signals, percent1, percent2)
    outputFile(metrics[0], metrics[1], metrics[2], output_path)
    removeDecompressFile(filename)

