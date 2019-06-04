import gzip
import re
import argparse
import os
import logging

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
	return file[:-3] #returns filename without .gz extension

def sortSignals(filename):
	"""
		Reads a wig or bedgraph file and sorts the signals by strongest to weakest.
	"""
	fileExt = getFileExtension(filename)
	#if file gziped we create a new decompressed file
	if fileExt == "gz":
		filename = decompressFile(filename)
		fileExt = getFileExtension(filename)

	if not (fileExt == "wig" or fileExt == "bedgraph"):
		raise Exception("Error : " + filename + " has to be .wig or .bedgraph")

	lines = [] #An array containing each lines of the file

	log.info("READING : " + os.path.basename(filename))
	with open(filename) as file:
		for line in file:
			parsed_line = line.split("\t")
			if not re.search("[a-zA-Z]", parsed_line[-1]): #checks if line is not a description line
				lines.append(parsed_line[-1])

	lines.sort(reverse = True)
	log.info("SORTING DONE !")
	return lines

def computeMetrics(sorted_signals):
	"""
		Calculates the sum of the signals, the top 10% and the top 5% bins.
	"""
	signalSum = 0
	sumTopBins_10 = 0
	sumTopBins_05 = 0
	
	#Creates all the values needed to calculate metrics
	nbBins = len(sorted_signals)
	nbTopBins_10 = int(round(nbBins * 0.1))
	nbTopBins_05 = int(round(nbBins * 0.05))

	for i in range(nbBins):
		signalSum += float(sorted_signals[i])

	for i in range(nbTopBins_10):
		if(i == nbTopBins_05):
			sumTopBins_05 = sumTopBins_10
		sumTopBins_10 += float(sorted_signals[i])
	return [signalSum, sumTopBins_10, sumTopBins_05]

def outputFile(signalSum, sumTopBins_10, sumTopBins_05, output_path):
	"""
		Outputs in a tab separated file the values of the signal sum, the sum of the top 10% bins,
		the sum of the top 5% bins, the ratio of the top 10% bins and the ratio of the top 5% bins.
	"""
	columns = "Signal sum\tTop 10% bins sum\tTop 5% bins sum\tRatio top 10% bins\tRatio top 5% bins"
	values = str(signalSum)+"\t"+str(sumTopBins_10)+"\t"+str(sumTopBins_05)+"\t"+str(sumTopBins_10/signalSum)+"\t"+str(sumTopBins_05/signalSum)

	complete_file_name = os.path.join(output_path, "signalToNoise.tsv")

	#Writes metrics in a new file
	log.info("WRITING METRICS IN : " + output_path)
	output = open(complete_file_name, "w+")
	output.write(columns+"\n")
	output.write(values)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Calculates the amount of noise in a signal track")
	parser.add_argument("-i", "--input_file", help="Path to input file to analyse", type=str, required=True)
	parser.add_argument("-o", "--output", help="Output directory", type=str, required=False, default=".")
	args = parser.parse_args()

	log = logging.getLogger("Download data")
	log.setLevel(logging.DEBUG)

	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)

	format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	ch.setFormatter(format)
	log.addHandler(ch)

	filename = args.input_file #Observed track
	output_path = args.output
	sorted_signals = sortSignals(filename)
	metrics = computeMetrics(sorted_signals)
	outputFile(metrics[0], metrics[1], metrics[2], output_path)
	log.info("RUN OVER")