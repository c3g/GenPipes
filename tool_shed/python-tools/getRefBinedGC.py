#!/usr/bin/python

### Mathieu Bourgey (2013/03/21)
### get GC and mappability count by bin

import os
import sys
import string
import getopt
import re
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from bedtools import IntervalFile, Interval


def getarg(argument):
	refF=""
	optli,arg = getopt.getopt(argument[1:],"s:r:m:o:h",['siz','ref','map','out','help'])
	if len(optli) < 4 :
		usage()
		sys.exit("Error : Missing argument(s)")
	for option, value in optli:
		if option in ("-s","--siz"):
			bsi=int(value)
		if option in ("-r","--ref"):
			refF=str(value)	
		if option in ("-m","--map"):
			mapp=str(value)	
		if option in ("-o","--output"):
			out=str(value)
		if option in ("-h","--help"):
			usage()
			sys.exit()
	if not os.path.exists(refF) :
		sys.exit("Error - reference file not found:\n"+refR)
	return  bsi, refF, mapp, out



def usage():
	print "USAGE : getRefBinedGC.py [option] "
	print "       -s :        bin size (bp)"
	print "       -r :        reference genome file "
	print "       -m :        reference mappability file "
	print "       -o :        output basename"
	print "       -h :        this help \n"


def main():
	print "\n---------------------------------------------------------------------------------"
	print "getRefBinedGC.py will generate a bed files with the GC of each bin from "
	print "all the reference genome "
	print "This program was written by Mathieu BOURGEY"
	print "For more information, contact: mbourgey@genomequebec.com"
	print "----------------------------------------------------------------------------------\n"
	siz, refF, mapF, outF = getarg(sys.argv)
	mapp= IntervalFile(mapF)
	out=open(outF,'w')
	out.write("Chrom\tStart\tEnd\tGCcontent\tUnMappContent\n")
	for seq_record in SeqIO.parse(refF, "fasta"):
		chroS=len(seq_record)
		chroID=str(seq_record.id)
		start=0
		end=start+siz-1
		while start <= chroS :
			if end > chroS:
				end=chroS
			i = Interval(chrom=chroID,start=start+1,end=end+1)
			refSeq=seq_record.seq[start:end]
			GCcont=int(GC(refSeq))
			h=mapp.all_hits(i)
			ct=0.0
			for i in range(0,len(h),1) :
				if h[i].start < start+1 :
					stp=start+1
				else :
					stp=h[i].start
				if h[i].end > end+1:
					enp=end+1
				else :
					enp=h[i].end
				ct=ct+int(enp-stp)
			MapCont=int((ct/float(siz))*10000)/100.00
			## write Output
			out.write(chroID +"\t"+str(start)+"\t"+str(end)+"\t"+str(GCcont)+"\t"+str(MapCont)+"\n")
			start=start+siz
			end=start+siz-1
	out.close()


main()

	    
	
	





