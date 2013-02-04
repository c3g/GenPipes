#!/usr/bin/python

import os
import sys
import string
import getopt
import re

class aliasPE:
	def __init__(self):
		self.R1=[]
		self.R2=[]
		self.runR1=""
		self.runR2=""
		self.mergeF=False
	
	def checkMultiF(self):
		##warning: run commands do not include output file !!
		doM=False
		if (len(self.R1) > 1) and (len(self.R2) > 1) :
			doM=True
			self.mergeF=True
			self.runR1="zcat "
			self.runR2="zcat "
			for j in range(0,len(self.R1),1) :
				self.runR1=self.runR1+self.R1[j]+" "
				self.runR2=self.runR2+self.R2[j]+" "
		return doM

	def merge(self,name,fi,outP,tmpout):
		##warning: run commands do not include output file !!
		if self.mergeF :
			self.runR1=self.runR1+"| gzip -c > "+outP+"/"+name+"_Merge_R1.fastq.gz\n"
			self.runR2=self.runR2+"| gzip -c > "+outP+"/"+name+"_Merge_R2.fastq.gz\n"
			tmpout.write(self.runR1)
			tmpout.write(self.runR2)
			fi.write(outP+"/"+name+"_Merge_R1.fastq.gz,"+name+"\n")
			fi.write(outP+"/"+name+"_Merge_R2.fastq.gz,"+name+"\n")
		else :
			fi.write(self.R1[0]+","+name+"\n")
			fi.write(self.R2[0]+","+name+"\n")
	

class aliasSI:
	def __init__(self):
		self.sing=[]
		self.runSI=""
		self.mergeF=False
		
	def checkMultiF(self):
		##warning: run commands do not include output file !!
		doM=False
		if len(self.sing) > 1 :
			doM=True
			self.runSI="zcat "
			for j in range(0,len(self.sing),1) :
				self.runSI=self.runSI+self.sing[j]+" "
		return doM

	def merge(self,name,fi,outP,tmpout):
		##warning: run commands do not include output file !!
		if self.mergeF :
			self.runSI=self.runSI+"| gzip -c > "+outP+"/"+name+"_Merge_SI.fastq.gz"
			tmpout.write(self.runSI)
			fi.write(outP+"/"+name+"_Merge_SI.fastq.gz,"+name+"\n")
		else :
			fi.write(self.sing[0]+","+name+"\n")


def getFastqListSI(x):
	fi=open(x,'r')
	li=fi.readline()
	dico={}
	while li != "" :
		ci=re.split('[,\n]',li)
		if not dico.has_key(ci[1]) :
			dico[ci[1]]=aliasSI()
		dico[ci[1]].sing.append(ci[0])
		li=fi.readline()
	fi.close()
	return dico
	

def getFastqListPE(x):
	fi=open(x,'r')
	li=fi.readline()
	dico={}
	while li != "" :
		ci=re.split('[,\n]',li)
		if not dico.has_key(ci[1]) :
			dico[ci[1]]=aliasPE()
		teR1=re.search("_R1",ci[0])
		if teR1 == None :
			dico[ci[1]].R2.append(ci[0])
		else :
			dico[ci[1]].R1.append(ci[0])
		li=fi.readline()
	fi.close()
	return dico
			
			
		
###Get argument form command line
def getarg(argument):
	optli,arg = getopt.getopt(argument[1:],"f:d:o:h",['file','output','design','help'])
	if len(optli) == 0 :
		usage()
		sys.exit("Error : No argument given")
	for option, value in optli:
		if option in ("-f","--file"):
			fil=str(value)
		if option in ("-o","--output"):
			if str(value)[-1] == "/" :
				out=str(value)[:-1]
			else :
				out=str(value)
		if option in ("-d","--design"):
			design=int(value)
		if option in ("-h","--help"):
			usage()
			sys.exit()
	return fil, out, design


###Usage information
def usage():
	print "\n-----------------------------------------------------------------------------"
	print "AliasFastqMerge.py checks the corresponding between fastyq file and the alias"
	print "If an alias has more than one fatstq file for each ends (paired-end) or for each"
	print "single read (single) this files will be merged, new file/alias file will be created"
	print "and the original filealais will be save as <NAME_FILE_ALIAS>.old"
	print "This program was written by Mathieu BOURGEY"
	print "For more information, contact: mbourgey@genomequebec.com"
	print "------------------------------------------------------------------------------\n"
	print "USAGE : AliasFastqMerge.py [option] "
	print "       -f :        alias file"
	print "       -o :        absolute path for the merge fatsq output folder"
	print "       -d :        design (1: Paired-end ; 0: Single)"
	print "       -h :        this help \n"



def main():
	fa, outP, des = getarg(sys.argv) ###Get argument form command line
	mc=[]
	if des == 0 :
		AFlist=getFastqListSI(fa)
	elif des == 1 :
		AFlist=getFastqListPE(fa)
	else :
		usage()
		sys.exit("Wrong design argument\n")
	for i in AFlist.keys() :
		mc.append(AFlist[i].checkMultiF())
	if sum(mc) >= 1 :
		os.system("cp "+fa+" "+fa+".old")
		fi=open(fa,'w')
		outtmp=open("runMergeFa.sh",'w')
		for i in AFlist.keys() :
			newAF=AFlist[i].merge(i,fi,outP,outtmp)
		fi.close()
		outtmp.close()
		print "0"
	else :
		print "1"
		

main()

