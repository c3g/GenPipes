#!/usr/bin/python

### Mathieu Bourgey (2013/08/15) - mbourgey@genomequebec.com
### format gtf obtained form cuffcompare (from cufflinks)


import os
import sys
import string
import getopt
import re

class transcript :  ## general transcript class that combine all information
	def __init__(self) :
		self.chro="" ## fill in tmpExon ini then addExon object method
		self.met=""  ## fill in tmpExon ini then addExon object method
		self.strand=""  ## fill in tmpExon ini then addExon object method
		self.geneID=""  ## fill in getTransTracking()
		self.nearRef=""  ## fill in tmpExon ini then addExon object method
		self.geneName=""  ## fill in getTransTracking()
		self.RefInfo=""  ## fill in getTransTracking()
		self.exSt=[]  ## fill in tmpExon ini then addExon object method
		self.exEn=[]## fill in tmpExon ini then addExon object method
		self.fpkm=[]  ## fill in getTransTracking()
		
	
	def addExon(self,tmpE): ## add infor from the tmpExon class
		if self.chro == "" :
			self.chro=tmpE.chro
			self.met=tmpE.met
			self.strand=tmpE.strand
			self.nearRef=tmpE.nearRef
		self.exSt.append(int(tmpE.st))
		self.exEn.append(int(tmpE.en))
	
	def getoutput(self,n,o): ## write transcript in a tsv output file (arg o)
		output=[]
		if len(self.exSt) >= 1 and len(self.exEn) >= 1 and len(self.exSt) == len(self.exEn) :
			stT=min(self.exSt)
			enT=max(self.exEn)
			exSL=",".join(map(str,self.exSt))
			exEL=",".join(map(str,self.exEn))
			output=[n,self.met,self.chro,str(stT),str(enT),self.strand,exSL,exEL,self.geneID,self.nearRef]
			output.extend(self.fpkm)
			output.append(self.RefInfo)
			o.write("\t".join(output)+"\n")


class tmpExon : ## exon class obtained when going over the column of the combined gtf file
	def __init__(self,c) :
		tmpC=re.sub(" ","",c[8])
		info=c[8].split(";")
		self.add=True
		self.nearRef=""
		for i in range(0,len(info),1):
			infoSpec=info[i].split("\"")
			classN=re.sub(" ","",infoSpec[0])
			if classN == "contained_in" :
				self.add=False
				break
			elif classN  == "nearest_ref" :
				self.nearRef=infoSpec[1]
			elif classN  == "gene_id" :
				self.geneID=infoSpec[1]
		if self.add :
			self.chro=c[0]
			self.met=c[1]
			self.st=int(c[3])
			self.en=int(c[4])
			self.strand=c[6]


def getarg(argument):
	c=""
	t=""
	s=""
	out=""
	optli,arg = getopt.getopt(argument[1:],"c:o:t:s:h",['comb','output','track','sample','help'])
	if len(optli) == 0 :
		usage()
		sys.exit("Error : No argument given")
	for option, value in optli:
		if option in ("-c","--comb"):
			c=str(value)
		if option in ("-t","--track"):
			t=str(value)
		if option in ("-s","--sample"):
			s=str(value)
		if option in ("-o","--output"):
			out=str(value)
		if option in ("-h","--help"):
			usage()
			sys.exit()
	if not os.path.exists(c) :
		sys.exit("Error - file not found:\n"+str(c))
	if not os.path.exists(t) :
		sys.exit("Error - file not found:\n"+str(t))
	if not os.path.exists(s) :
		sys.exit("Error - file not found:\n"+str(s))
	return c, t, s, out

def usage():
	print "\n---------------------------------------------------------------------------------"
	print "formatDenovoCombinedGTF.py format gtf obtained form cuffcompare (from cufflinks)"
	print "This program was written by Mathieu BOURGEY"
	print "For more information, contact: mbourgey@genomequebec.com"
	print "----------------------------------------------------------------------------------\n"
	print "USAGE : GetRecurenceSV.py [option] "
	print "       -c :        combined gtf"
	print "       -t :        tracking file"
	print "       -s :        sample list"
	print "       -o :        output gtf"
	print "       -h :        this help \n"

def getCode():
	dico={}
	dico["="]="Complete match of intron chain"
	dico["c"]="Contained"
	dico["j"]="Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript"
	dico["e"]="Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment"
	dico["i"]="A transfrag falling entirely within a reference intron"
	dico["o"]="Generic exonic overlap with a reference transcript"
	dico["p"]="Possible polymerase run-on fragment (within 2Kbases of a reference transcript)"
	dico["r"]="Repeat. Currently determined by looking at the soft-masked reference sequence and applied to transcripts where at least 50% of the bases are lower case"
	dico["u"]="Unknown, intergenic transcript"
	dico["x"]="Exonic overlap with reference on the opposite strand"
	dico["s"]="An intron of the transfrag overlaps a reference intron on the opposite strand (likely due to read mapping errors)"
	dico["."]="(.tracking file only, indicates multiple classifications)"
	return dico
	

def getSampleOrder(x):
	f=open(x,'r')
	l=f.readline()
	sample=[]
	while l != "" :
		co=l.split()
		sample.append("FPKM_"+co[0])
		l=f.readline()
	order="\t".join(sample)
	f.close()
        return order

def getTransTracking(x) :
	f=open(x,'r')
	l=f.readline()
	codeC=getCode()
	trans={}
	while l != "" :
		c=l.split("\t")
		try :
			c[-1]=re.sub("\n","",c[-1]) 
		except :
			stop("tracking file not in unix text file format try dos2Unix\n")
		trans[c[0]]=transcript()
		trans[c[0]].geneID=c[1]
		ref=c[2].split("|")
		if len(ref) > 1 :
			trans[c[0]].geneName=ref[0]
		else :
			trans[c[0]].geneName=c[2]
		trans[c[0]].RefInfo=codeC[c[3]]
		for i in range(4,len(c),1) :
			if c[i] != "-" :
				si=c[i].split("|")
				trans[c[0]].fpkm.append(si[3])
			else :
				trans[c[0]].fpkm.append(c[i])
		
		l=f.readline()
	f.close()
	return trans

def matchTransComb(x,trans) :
	f=open(x,'r')
	l=f.readline()
	exon={}
	while l != "" :
		c=l.split("\t")
		try :
			c[-1]=re.sub("\n","",c[-1]) 
		except :
			stop("combined gtf file not in unix text file format try dos2Unix\n")
		name=l.split("\"")
		exon[name[3]]=tmpExon(c)
		if exon[name[3]].add :
			trans[name[3]].addExon(exon[name[3]])
		l=f.readline()
	f.close()
	return trans


      
def main():
	comb, track, samp, out =  getarg(sys.argv)
	sampOrd=getSampleOrder(samp)  ## Get sample order for header 
	trans=getTransTracking(track)   ## Get transcript basic info
	trans=matchTransComb(comb,trans)  ## add exon info to trnascripts
	transID=trans.keys()
	transID.sort()
	of=open(out,'w')
	header=["transcript_ID","method","chromosome","start","end","strand","exon_starts","exon_ends","gene_ID","nearest_ref",sampOrd,"reference_match"]
	of.write("\t".join(header) + "\n")
	for i in transID:
		trans[i].getoutput(i,of) ##output transcipts
	of.close()

main()

