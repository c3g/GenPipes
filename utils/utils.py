#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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
import argparse
import glob
import logging
import os
import re
import string
import sys
import shutil
import subprocess
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
        if x[-1].upper() in dico:
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
        print(i +"is now removed\n")


def time_to_datetime(time):
    """
        From transform time str to delta time:
        Acceptable time formats include "minutes", "minutes:seconds",
        "hours:minutes:seconds", "days-hours",
        "days-hours:minutes" and "days-hours:minutes:seconds"
        In fact it will also work for pbs/torque.
        Note that is strip from the string everithing that is not the time
    :param time: sting containing a slurm "--time" format substring
    :return: timedelta object
    """

    time = re.search("([0-9]+-)?[0-9]+:[0-9]+(:[0-9]+)?", time).group()

    if '-' in time:
        days, rest = time.split('-')
        rest = rest.split(':')[0]
    else:
        rest = time.split(':')
        days = 0

    hours = int(rest[0]) + int(days * 24)
    minutes = int(rest[1])

    if len(rest) > 2:
        sec = int(rest[2])
    else:
        sec = 0

    return datetime.timedelta(hours=hours, minutes=minutes, seconds=sec)


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


def container_wrapper_argparse(argv):
    """
        This method start the pipeline inside de system put in place by the genpipe_in_a_container
        cvmfs/container wrapper.

    Args:
        argv: Any valid Genpipes option plus --wrap.
    Returns: Return code of the subprocess

    """

    help = False
    if '-h' in argv:
        try:
            argv.remove('-h')
        except ValueError:
            pass
        help = True

    script_dir_current = os.path.dirname(os.path.realpath(__file__))
    parser = argparse.ArgumentParser(conflict_handler='resolve')
    parser.add_argument('--wrap', type=str, help="Path to the genpipe cvmfs wrapper script.\n"
                                                 "Default is genpipes/ressources/container/bin/container_wrapper.sh",
                        nargs='?')
    args, argv = parser.parse_known_args(argv)
    if args.wrap is None:
        genpipes_home = '/'.join(script_dir_current.split('/')[:-1])
        default_wrapper = '{}/resources/container/bin/container_wrapper.sh'.format(genpipes_home)
        args.wrap = default_wrapper

    wrap_option = ['--container', 'wrapper', args.wrap]

    if (args.wrap and 'batch' in argv) and '--no-json' not in argv:
        parser.error("--wrap and -j batch requires --no-json")

    sys.stderr.write('wrapping\n')
    # call in the wrapper
    if help:
        argv.append('--help')
    sys.stderr.write('{} {} {}'.format(args.wrap, ' '.join(argv), ' '.join(wrap_option)))
    return subprocess.call([args.wrap] + argv + wrap_option)

