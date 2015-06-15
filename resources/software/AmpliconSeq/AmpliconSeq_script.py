#!/usr/bin/env python
# -*- coding: utf-8 -*- 


################################################################################
### Several scripts for the Amplicon-Seq pipeline including outputs creation.
################################################################################

import sys, getopt
import os, re

#m: map_build, krona

def map_build(samples):

	out_map = open('map.txt','w')
	out_map.write("#SampleID\tInputFileName\tDescription\n")
	
	for sample in samples.split(','):
		out_map.write(sample+"\t"+sample+".fastq\t"+sample+"\n")

def krona(table_tax):

	table = open(table_tax,"r")	
	lines = table.readlines()
	
	word = lines[1].split()
	sample_name = word[2:-1]
	file_name=[]
	
	sample_number = 0		
		
	while sample_number < len(sample_name):
	
		out_krona = open("alpha_diversity/krona_chart/"+sample_name[sample_number]+".txt","w")
		file_name.append("alpha_diversity/krona_chart/"+sample_name[sample_number]+".txt")
		
		i=2
		while i < len(lines):
		
			word = lines[i].split()	
			out_krona.write(word[sample_number+1]+'\t'+'\t'.join(word[len(sample_name)+1:])+"\n")
			i+=1
		
		sample_number += 1	
		out_krona.close()
	
	table.close()
	

def main(argv):
	inputfile = ''
	outputfile = ''
	
	opts, args = getopt.getopt(argv,"hi:o:m:")
	if len(opts) == 0 and len(args) == 0:
		print 'test.py -i <inputfile> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'test.py -i <inputfile> -o <outputfile>'
			sys.exit()
		elif opt in ("-m"):
			if arg == "map_build":
			
				input_map = False
				for optM, argM in opts:					
					if optM in ("-i"):						
						input_map = True
						break
					
				if input_map:
					print("Map file building ...")
					map_build(argM)
				else:
					print("Need input")
					
			elif arg == "krona":
			
				input_krona = False
				for optM, argM in opts:					
					if optM in ("-i"):						
						input_krona = True
						break
					
				if input_krona:
					print("Krona chart building ...")
					krona(argM)
				else:
					print("Need input")


			else :
					print("Argument availaible for -m: map_build, krona")

					
		elif opt in ("-o", "--ofile"):
			outputfile = arg

if __name__ == "__main__":
	main(sys.argv[1:])
	