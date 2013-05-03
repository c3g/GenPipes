#!/usr/bin/python

### Mathieu Bourgey (2012/01/04) - mbourgey@genomequebec.com

import os
import sys
import string
import getopt
import re


def getarg(argument):
	fa=""
	leng=""
	thr=1000
	dire="M"
	out=""
	optli,arg = getopt.getopt(argument[1:],"f:l:t:d:o:h",['fasta','length','threshold','direction','output','help'])
	if len(optli) == 0 :
		usage()
		sys.exit("Error : No argument given")
	for option, value in optli:
		if option in ("-f","--fasta"):
			if os.path.exists(str(value)) :
				fa=str(value)
			else :
				sys.exit("Error - fasta file not found:\n"+str(value))
		if option in ("-l","--length"):
			if os.path.exists(str(value)) :
				leng=str(value)
			else :
				sys.exit("Error - length file not found:\n"+str(value))
		if option in ("-t","--threshold"):
			thr=int(value)	
		if option in ("-d","--direction"):
			dire=str(value)
		if option in ("-o","--output"):
			out=str(value)
		if option in ("-h","--help"):
			usage()
			sys.exit()
	if out == "" :
		out=fa+".sizeSel.fa"
	if dire.upper() != "M" :
		dire="M"
	setPar(fa, leng, thr, dire, out)
	return fa, leng, thr, dire, out

def usage():
	print "\n-----------------------------------------------------------------------------------------------------"
	print "fasta_Length_select.py select a subset of fasta sequence based on a given sequence size ditribution"
	print "This program was written by Mathieu BOURGEY"
	print "For more information, contact: mbourgey@genomequebec.com"
	print "------------------------------------------------------------------------------------------------------\n"
	print "USAGE : fasta_Length_select.py [options] "
	print "mandatory options:"
	print "       -f :        fasta file"
	print "       -l :        length file"
	print "other options:"
	print "       -t :        threshold selection size in bp (default : 1000)"
	print "       -d :        direction of selection (m/M) : m - to lower ; M - to upper (default: M)"
	print "       -o :        output file (default: <fasta file>.sizeSel.fa"
	print "       -h :        this help \n"
	
def setPar(f, l, t, d, o):
	dico={}
	dico["m"]="to lower"
	dico["M"]="to upper"
	print "parameters:"
	print "fasta file  : "+f
	print "length file : "+l
	print "output file : "+o
	print "threshold   : "+str(t)
	print "direction   : "+dico[d]
	
## get list of the name of fatsa sequence to keep
def getList(l,t,d):
	print "Generate inclusion list... "
	listFaUp=[]
	listFaDown=[]
	fi=open(l,'r')
	li=fi.readline()
	while li != "" :
		ci=li.split()
		if int(ci[0]) >= t :
			listFaUp.append(">"+ci[1])
		if int(ci[0]) <= t : 
			listFaDown.append(">"+ci[1])
		li=fi.readline()
	fi.close()
	print "to upper : "+str(len(listFaUp))+" sequences"
	print "to lower : "+str(len(listFaDown))+" sequences"
	print "Done\n"
	if d == "m" :
		return listFaDown
	else :
		return listFaUp
	

def main():
	ff, lf, th, di, outf = getarg(sys.argv)
	fkeep=getList(lf,th,di)
	fi=open(ff,'r')
	out=open(outf,'w')
	li=fi.readline()
	ow=False
	print "Generate output ... "
	while li != "" :
		if li[0] == ">" :
			ow=False
			ci=li.split()
			if ci[0] in fkeep :
				ow=True
		if ow :
			out.write(li)
		li=fi.readline()
	fi.close()
	out.close()
	print "Done\n"
	
main()