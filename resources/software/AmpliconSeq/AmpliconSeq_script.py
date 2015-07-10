#!/usr/bin/env python
# -*- coding: utf-8 -*- 


################################################################################
### Several scripts for the Amplicon-Seq pipeline including outputs creation.
################################################################################

import sys, getopt
import os, re

#m: map_build, krona

def map_build(samples):

	"""
	Make a general map file if not available.
	"""

	out_map = open('map.txt','w')
	out_map.write("#SampleID\tInputFileName\tDescription\n")
	
	for sample in samples.split(','):
		out_map.write(sample+"\t"+sample+".fastq\t"+sample+"\n")
		
	out_map.close()
	
def map_per_sample(sample,output):

	"""
	Make a basic map file for each sample.
	"""

	out_map = open(output,'w')
	out_map.write("#SampleID\t"+sample+"\n")
	out_map.write(sample+"\t"+sample+"\n")
		
	out_map.close()
	

def krona(table_tax):

	"""
	Create file for Krona chart.
	"""

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

def uchime_stat(uchime_log, flash_log, sample):

	"""
	Uchime statistics for report.
	"""

	log_merged = open(flash_log,"r")
	log_chimer = open(uchime_log,"r")

	filter_stat=[]
	filter_stat.append([int(i.split()[3]) for i in log_merged if re.search("Combined pairs",i)][0])
	filter_stat.append([int(float(i.split()[1])) for i in log_chimer if re.search(sample,i)][0])
	filter_stat.append(round(float(filter_stat[1])/float(filter_stat[0])*100,2))

	log_merged.close()
	log_chimer.close()
	
	#print filter_stat
	print "\t".join(str(x) for x in filter_stat)
	
def sample_rarefaction(alpha_stat_f, alpha_out_f, sample_name):

	"""
	Create rarefaction files for differents metrics 
	(observed species, chao1 and shannon) for each sample.
	"""
	
	alpha_stat = open(alpha_stat_f,'r')
	alpha_out = open(alpha_out_f,'w')
	
	lines = alpha_stat.readlines()
	
	word = re.split("[\r\t\n]",lines[0])
	index_sample = word.index(sample_name)
	
	i=0
	while i<len(lines):
		word = re.split("[\r\t\n]",lines[i])
		alpha_out.write('\t'.join(word[0:3])+'\t'+word[index_sample]+'\n')
		i+=1
	
	alpha_stat.close()
	alpha_out.close()
	
def single_rarefaction(alpha_stat_f, alpha_out_f, rarefaction_threshold):

	"""
	Rerefying all sample at the same level.
	"""

	alpha_stat = open(alpha_stat_f,'r')
	alpha_out = open(alpha_out_f,'w')
	
	lines = alpha_stat.readlines()
	rarefaction_done = False
	
	i=0
	
	while i<len(lines):
		word = re.split("[\r\t\n]",lines[i])
		
		if rarefaction_done == False:
			if word[0] == """alpha_rarefaction_{}_2.txt""".format(int(rarefaction_threshold)):			
				rarefaction_done = True
				
			alpha_out.write(lines[i])
			i+=1
			
		else:
			break
	
	alpha_stat.close()
	alpha_out.close()

def plot_heatmap(table_f, rep_out_f, taxon_lvl):

	"""
	Create R script for heatmap plot.
	"""

	if int(taxon_lvl) == 0:
		name_tax = "Domain"
		colnames_for_R = 'colnames(taxmat) <- c("Domain")'
	
	if int(taxon_lvl) == 1:
		name_tax = "Phylum"
		colnames_for_R = 'colnames(taxmat) <- c("Domain", "Phylum")'
	
	if int(taxon_lvl) == 2:
		name_tax = "Class"
		colnames_for_R = 'colnames(taxmat) <- c("Domain", "Phylum", "Class")'
		
	if int(taxon_lvl) == 3:
		name_tax = "Order"
		colnames_for_R = 'colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order")'
		
	if int(taxon_lvl) == 4:
		name_tax = "Family"
		colnames_for_R = 'colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family")'
		
	if int(taxon_lvl) == 5:
		name_tax = "Genus"
		colnames_for_R = 'colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")'		
	
	if int(taxon_lvl) == 6:
		name_tax = "Species"
		colnames_for_R = 'colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")'	
				
		
	known_phylum=['Acidobacteria', 'Actinobacteria', 'Aquificae', 'Armatimonadetes', 'Bacteroidetes', 'Caldiserica', 'Chlamydiae', 'Chlorobi', 'Chloroflexi', 'Chrysiogenetes', 'Cyanobacteria', 'Deferribacteres', 'Deinococcus-Thermus', 'Dictyoglomi', 'Elusimicrobia', 'Fibrobacteres', 'Firmicutes', 'Fusobacteria', 'Gemmatimonadetes', 'Lentisphaerae', 'Nitrospira', 'Planctomycetes', 'Proteobacteria', 'Spirochaetes', 'Synergistetes', 'Tenericutes', 'Thermodesulfobacteria', 'Thermomicrobia', 'Thermotogae', 'Verrucomicrobia','Unclassified']
		
	#1st step: Create the OTU matrix
	
	table = open(table_f,"r")
	
	lines = table.readlines()
	
	row_names = []
	data = []
	OTU_tax_final = []
	
	
	header_sample = lines[1].split()
	col_names = header_sample[2:]	#Name of sample
	
	num_sampe = len(col_names)		#Number of sample
	
	i = 2
	while i < len(lines):
	
		parse = lines[i].split()
		
		for taxon in parse[0].split(';'):
			if taxon == 'Other' and OTU_tax_final[len(OTU_tax_final)-1]!='k__Bacteria':
				OTU_tax_final.append(str(OTU_tax_final[len(OTU_tax_final)-1])+'_'+taxon)
			else:
				OTU_tax_final.append(taxon)
		
		j = 1
		while j<len(parse):
			data.append(parse[j])
			j+=1	
		
		row_names.append("Otu{}".format(i-1))
											
		i+=1
	
	table.close()
	
	out_subprocess = open(os.path.join(rep_out_f,"OTU_data.txt"),"w")
	for OTU_data in data:
		out_subprocess.write(OTU_data+"\n")
	out_subprocess.close()
	
	out_subprocess = open(os.path.join(rep_out_f,"OTU_name.txt"),"w")
	for OTU_name in row_names:
		out_subprocess.write(OTU_name+"\n")
	out_subprocess.close()
	
	
	#2nd step: Create the taxonomy matrix
	
	
	out_subprocess = open(os.path.join(rep_out_f,"OTU_tax_final.txt"),"w")
	for taxon in OTU_tax_final:
		out_subprocess.write(taxon+"\n")
	out_subprocess.close()
	
	
	#######################################################
	
	
	cmd_matrix = "otumat = matrix(vec_otu_data, nrow={}, ncol={}, byrow=TRUE)".format(len(row_names), len(col_names))
	cmd_names = "dimnames(otumat) = list(vec_otu_name, c({}))".format(str(col_names)[1:len(str(col_names))-1])
	
	out_to_R = open(os.path.join(rep_out_f,"OTU_%s_to_R.R" % (name_tax)),"w")
	
	out_to_R.write('#!/usr/bin/Rscript\n')
	out_to_R.write('library("pheatmap")\n')
	out_to_R.write('library("RColorBrewer")\n')
	out_to_R.write('## if not installed, quickly add it as follows:\n')
	out_to_R.write('#source("http://bioconductor.org/biocLite.R")\n')
	out_to_R.write('#biocLite(c("RColorBrewer", "pheatmap"))\n\n')
		
	out_to_R.write("#Data load:\n\n")
	
	out_to_R.write('otu_data <- read.table("'+os.path.join(rep_out_f,'OTU_data.txt')+'")'+"\n")
	out_to_R.write('vec_otu_data<-as.vector(as.matrix(otu_data))'+"\n")
	out_to_R.write(cmd_matrix+"\n")
	out_to_R.write('otu_name <- read.table("'+os.path.join(rep_out_f,'OTU_name.txt')+'")'+"\n")
	out_to_R.write('vec_otu_name<-as.vector(as.matrix(otu_name))'+"\n")
	out_to_R.write(cmd_names+"\n")
	out_to_R.write('tax <- read.table("'+os.path.join(rep_out_f,'OTU_tax_final.txt')+'")'+"\n")
	out_to_R.write('vec_tax<-as.vector(as.matrix(tax))'+"\n")
	out_to_R.write('taxmat = matrix(vec_tax, nrow = nrow(otumat), ncol = %s,byrow=TRUE)' % (int(taxon_lvl)+1)+"\n")
	out_to_R.write('rownames(taxmat) <- rownames(otumat)'+"\n\n")
	out_to_R.write(colnames_for_R+"\n")
	
	out_to_R.write('#Heatmap:\n\n')
	out_to_R.write('png("'+os.path.join(rep_out_f,'otu_heatmap.png')+'",width = 800, height = 800, units = "px")\n')
	out_to_R.write('pheatmap(otumat,labels_row=taxmat[,"%s"])'% name_tax+'\n\n')
	out_to_R.write('dev.off()\n\n')
	
	out_to_R.write('write.table(otumat, "'+os.path.join(rep_out_f,'otumat.tsv')+'", sep="\t")\n')
	out_to_R.write('write.table(taxmat, "'+os.path.join(rep_out_f,'taxmat.tsv')+'", sep="\t")\n')
	
	out_to_R.close()
	

def main(argv):
	inputfile = ''
	outputfile = ''
	
	opts, args = getopt.getopt(argv,"hi:o:m:j:s:")
	if len(opts) == 0 and len(args) == 0:
		print("AmpliconSeq_script.py -m <arg>")
		print("Argument availaible for -m: map_build, map_per_sample, krona, uchime, sample_rarefaction, single_rarefaction")
		sys.exit(2)
		
	for opt, arg in opts:
		if opt == '-h':
			print("AmpliconSeq_script.py -m <arg>")
			print("Argument availaible for -m: map_build, map_per_sample, krona, uchime, sample_rarefaction, single_rarefaction")
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
					print("AmpliconSeq_script.py -m map_build -i <samples>")
					
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
					print("AmpliconSeq_script.py -m krona -i <krona_file>")
					
			elif arg == "uchime":
			
				input_uchime = False
				for optM, argM in opts:					
					if optM in ("-i"):	#uchime log
						uchime_log = argM
					
						for optM2, argM2 in opts:	#flash log				
							if optM2 in ("-j"):		
								flash_log = argM2
						
								for optM3, argM3 in opts:	#sample name				
									if optM3 in ("-s"):	
										sample = argM3.replace('_','.')
														
										input_uchime = True
										break
					
				if input_uchime:
					uchime_stat(uchime_log, flash_log, sample)
				else:
					print("AmpliconSeq_script.py -m uchime -i <uchime_log> -j <flash_log> -s <samples>")

			elif arg == "sample_rarefaction":
			
				input_sample_rarefaction = False
				for optM, argM in opts:					
					if optM in ("-i"):	#stat metrics
						alpha_stat_f = argM
					
						for optM2, argM2 in opts:	#output			
							if optM2 in ("-j"):		
								alpha_out_f = argM2
						
								for optM3, argM3 in opts:	#sample name				
									if optM3 in ("-s"):	
										sample = argM3.replace('_','.')
														
										input_sample_rarefaction = True
										break
					
				if input_sample_rarefaction:
					sample_rarefaction(alpha_stat_f, alpha_out_f, sample)
				else:
					print("AmpliconSeq_script.py -m sample_rarefaction -i <stat_metrics> -j <output_file> -s <samples>")
			
			elif arg == "map_per_sample":
			
				input_map_per_sample = False
				
				for optM2, argM2 in opts:	#flash log				
					if optM2 in ("-j"):		
						output = argM2
				
						for optM3, argM3 in opts:	#sample name				
							if optM3 in ("-s"):	
								sample = argM3.replace('_','.')
												
								input_map_per_sample = True
								break
					
				if input_map_per_sample:
					map_per_sample(sample,output)
				else:
					print("AmpliconSeq_script.py -m map_per_sample -s <samples> -j <output_file>")
					
			elif arg == "single_rarefaction":
			
				input_single_rarefaction = False
				for optM, argM in opts:					
					if optM in ("-i"):	#stat metrics
						alpha_stat_f = argM
					
						for optM2, argM2 in opts:	#output			
							if optM2 in ("-j"):		
								alpha_out_f = argM2
						
								for optM3, argM3 in opts:	#rarefaction threshold				
									if optM3 in ("-s"):	
										rarefaction_threshold = argM3
														
										input_single_rarefaction = True
										break
					
				if input_single_rarefaction:
					single_rarefaction(alpha_stat_f, alpha_out_f, rarefaction_threshold)
				else:
					print("AmpliconSeq_script.py -m single_rarefaction -i <stat_metrics> -j <output_file> -s <rarefaction_threshold>")		
					
			elif arg == "plot_heatmap":
			
				input_plot_heatmap = False
				for optM, argM in opts:					
					if optM in ("-i"):	#table
						table_f = argM
					
						for optM2, argM2 in opts:	#output			
							if optM2 in ("-j"):		
								rep_out_f = argM2
						
								for optM3, argM3 in opts:	#taxon level			
									if optM3 in ("-s"):	
										taxon_lvl = argM3
														
										input_plot_heatmap = True
										break
					
				if input_plot_heatmap:
					plot_heatmap(table_f, rep_out_f, taxon_lvl)
				else:
					print("AmpliconSeq_script.py -m plot_heatmap -i <table.txt> -j <output_directory> -s <taxon_lvl>")	

			else :
					print("Argument availaible for -m: map_build, map_per_sample, krona, uchime, sample_rarefaction, single_rarefaction, plot_heatmap")

					
		elif opt in ("-o", "--ofile"):
			outputfile = arg

if __name__ == "__main__":
	main(sys.argv[1:])
	