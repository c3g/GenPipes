#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import ConfigParser
import glob
import logging
import os
import re
import string
import sys
import shutil
import datetime


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


def slurm_time_to_datetime(time):
    """
        From slurm doc:
        Acceptable time formats include "minutes", "minutes:seconds",
        "hours:minutes:seconds", "days-hours",
        "days-hours:minutes" and "days-hours:minutes:seconds"
    :param time: sting from slurm sbatch --time option
    :return: timedelta object
    """
    time = time.lstrip()
    time = time.lstrip('--time=')

    colon = time.count(':')
    dash = time.count('-')
    split_t = time.split(':')
    if len(split_t) == 1:
        split_d = split_t[0].split('-')
        if len(split_d) == 2:
            return datetime.timedelta(days=int(split_d[0]), hours=int(split_d[1]))
        return datetime.timedelta(minutes=int(split_t[0]))
    elif len(split_t) == 2:
        split_d = split_t[0].split('-')
        if len(split_d) == 2:
            return datetime.timedelta(days=int(split_d[0]),
                                      hours=int(split_d[1]), minutes=int(split_t[1]))
        return datetime.timedelta(hours=int(split_t[0]), minutes=int(split_t[1]))
    elif len(split_t) == 3:
        split_d = split_t[0].split('-')
        if len(split_d) == 2:
            return datetime.timedelta(days=int(split_d[0]),
                                      hours=int(split_d[1]),
                                      minutes=int(split_t[1]), seconds=int(split_t[2]))
        return datetime.timedelta(hours=int(split_t[0]),
                                  minutes=int(split_t[1]), seconds=int(split_t[2]))


def expandvars(path, skip_escaped=False):
    """Expand environment variables of form $var and ${var}.
       If parameter 'skip_escaped' is True, all escaped variable references
       (i.e. preceded by backslashes) are skipped.
       Unknown variables are set to '' like in a bash shell.
    """
    def replace_var(m):
        return os.environ.get(m.group(2) or m.group(1), '')
    reVar = (r'(?<!\\)' if skip_escaped else '') + r'\$(\w+|\{([^}]*)\})'
    return re.sub(reVar, replace_var, path)