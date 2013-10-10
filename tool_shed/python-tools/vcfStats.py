

### Mathieu Bourgey (2013/10/03)
### vcfStats

import os
import sys
import string
import getopt
import re
import vcf
import numpy as np
import matplotlib.pyplot as plt
#from Bio import SeqIO
#from Bio.Alphabet import IUPAC
#from Bio.SeqUtils import GC
#from bedtools import IntervalFile, Interval



def getarg(argument):
	vcfF=""
	dico="null"
	optli,arg = getopt.getopt(argument[1:],"v:d:o:h",['vcf','dico','output','help'])
	for option, value in optli:
		if option in ("-v","--vcf"):
			vcfF=str(value)
		if option in ("-d","--dico"):
			dico=str(value)	
		if option in ("-o","--output"):
			out=str(value)
		if option in ("-h","--help"):
			usage()
			sys.exit()
	if not os.path.exists(vcfF) :
		sys.exit("Error - vcf file not found:\n"+vcfF)
	return vcfF, dico, out



def usage():
	print "\n---------------------------------------------------------------------------------"
	print "vcfStats.py will generate SNV statistics from a vcf file "
	print "if no dico is given vcfStats will proceed for every chromosome "
	print "This program was written by Mathieu BOURGEY"
	print "For more information, contact: mbourgey@genomequebec.com"
	print "----------------------------------------------------------------------------------\n"
	print "USAGE : getRefBinedGC.py [option] "
	print "       -v :        vcf file"
	print "       -d :        dictionary file (optional)"
	print "       -o :        output basename"
	print "       -h :        this help \n"

def getRegion(x,r):
	d=[]
	l=x.readline()
	while l != "" :
		c=l.split()
		if c[0] == "@SQ" :
			chrI=c[1].split(":")
			if chrI[1] in r:
				d.append(chrI[1])
			else :
				print chrI[1] +" contigs not found: excluded from the analysis\n"
		l=x.readline()
	return d

#def generateGraphs(region,sample,out):
	####bar plot exemple to arrange 
	#fig = plt.figure()
	#ax = fig.add_subplot(111)
	### necessary variables
	#ind = np.arange(len(region.keys()))                # the x locations for the groups
	#width = 0.35                      # the width of the bars
	### the bars
	#rect=[]
	#values=[]
	#rk=region.keys()
	#maxV=0
	#for i in range(0,len(sample),1):
		#for j in range(0,len(rk),1):
			#if i == 0 :
				#values.append([])
			#if  int(region[j][i]) > maxV:
				#maxV=int(region[j][i])
			#values[i].append(int(region[j][i]))
	#colors = cm.rainbow(np.linspace(0, 1, len(sample)))
	#for i in range(0,len(sample),1) :
		#rect.append(ax.bar(ind,values[i] , width,color=colors[i]))
	## axes and labels
	#ax.set_xlim(-width,len(ind)+width)
	#ax.set_ylim(0,maxV*0.2)
	#ax.set_ylabel('Mutation rate')
	#ax.set_title('Mutation rate by chromosome and by sample')
	#xTickMarks = [i for i in rk]
	#ax.set_xticks(ind+width)
	#xtickNames = ax.set_xticklabels(xTickMarks)
	#plt.setp(xtickNames, rotation=45, fontsize=10)
	
	### add a legend
	#ax.legend( tuple(rect), tuple(sample) )
	#plt.savefig(out+"_SampleMutationRate.png", format='png')
	#plt.savefig(out+"_SampleMutationRate.pdf", format='pdf')

def main():
	vcfF, dicF, outB = getarg(sys.argv)
	vcf_reader = vcf.Reader(open(vcfF,'r'))
	if dicF != "null" and os.path.exists(dicF) :
		contigList=getRegion(open(dicF,'r'),vcf_reader.contigs.keys())
	else :
		print "No valid dictionary given: working the full contig list\n"
		contigList=vcf_reader.contigs.keys()
	if len(contigList) == 0 :
		print "Contigs in the given dictionary do not match any of those in the vsf: working the full contig list\n"
		contigList=vcf_reader.contigs.keys()
	region={}
	samples=vcf_reader.samples
	for record in vcf_reader:
		if len(record.REF) == 1 : ## only look at SNPs
			chr=record.CHROM
			if chr in contigList:
				if region.has_key(chr) :
					for i in samples:
						if record.genotype(i)['GT'] != "0/0" :
							region[chr][i]=region[chr][i]+1
				else :
					region[chr]={}
					for i in samples:
						if record.genotype(i)['GT'] != "0/0" :
							region[chr][i]=1
						else :
							region[chr][i]=0
	out=open(outB,'w')
	out.write("Contig\t"+"\t".join(region[region.keys()[0]].keys())+"\n")
	for i in region.keys() :
		clen=vcf_reader.contigs[i][1]
		for j in region[i].keys():
			region[i][j]=str(int(round(clen/region[i][j])))
		out.write(i+"\t"+"\t".join(region[i].values())+"\n")
	out.close()
	### generate the graph (png and pdf) for the report
	#generateGraphs(region,samples,outB)

main()

