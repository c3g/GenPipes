import os
import logging
import argparse


"""
Creates and input info file for chromimpute from a folder containing bigwig files.
"""

def createInputInfo(path_to_folder, output):
	filename = os.path.join(output, "inputinfo.txt")
	output = open(filename, "w+")
	log.info("========WRITING========")
	for filename in os.listdir(path_to_folder):
		parsed_file = filename.split(".")
		sample = parsed_file[2]
		mark = parsed_file[-4]
		output.write(sample+"\t"+mark+"\t"+filename[5:]+"\n")
	output.close()
	log.info("=========DONE==========")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Creates the input info file for ChromImpute")
	parser.add_argument("-p", "--path", help="Path to the folder containing the dataset to impute", type=str, required=True)
	parser.add_argument("-o", "--output", help="Path to the output directory", type=str, required=False, default=".")
	args = parser.parse_args()

	log = logging.getLogger("Download data")
	log.setLevel(logging.DEBUG)

	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)

	format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	ch.setFormatter(format)
	log.addHandler(ch)

	#TESTING
	path_to_folder = args.path
	output = args.output
	createInputInfo(path_to_folder, output)
