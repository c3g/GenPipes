#!/usr/bin/env python

# Python Standard Modules
import ConfigParser
import glob
import logging
import os
import re
import string
import sys
import shutil


# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

def number_symbol_converter(x):
    dico={}
    dico["K"]="000"
    dico["M"]="000000"
    dico["G"]="000000000"
    dico["T"]="000000000000"
    dico["P"]="000000000000000"
    dico["E"]="000000000000000000"
    dico["Z"]="000000000000000000000"
    dico["Y"]="000000000000000000000000"
    try:
        int(x[:-1])+1
        if dico.has_key(x[-1].upper()):
            return x[:-1]+dico[x[-1].upper()]
        elif int(x[-1])+1 :
            return x
        else :
            raise Exception("Number abbreviation \"" + x + "\" is not recognized finishing with a number symbol (k, M, G, T, P, E, Z or Y")
    except: 
        raise Exception("Number abbreviation \"" + x + "\" is not a number abbreviation")


def fpkm_correlation_matrix(cuffnorm_file, output_file):

    return Job(
        [cuffnorm_file],
        [output_file],
        [
            ['fpkm_correlation_matrix', 'module_R'],
        ],
        command="""\
R --no-save --no-restore <<-EOF
dataFile=read.table("{cuffnorm_file}",header=T,check.names=F)
fpkm=cbind(dataFile[,2:ncol(dataFile)])
corTable=cor(log2(fpkm+0.1))
corTableOut=rbind(c('Vs.',colnames(corTable)),cbind(rownames(corTable),round(corTable,3)))
write.table(corTableOut,file="{output_file}",col.names=F,row.names=F,sep="\t",quote=F)
print("done.")

EOF""".format(
        cuffnorm_file=cuffnorm_file,
        output_file=output_file,
    ))


def cleanFiles(x):
    ##x must be a list of files path 
    for i in x:
        shutil.rmtree(i)
        print i +"is now removed\n"

