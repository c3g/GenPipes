#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

# Python module where to install libs
PYTHON_VERSION=3.11.1
PYTHON_MODULE=mugqic/python/$PYTHON_VERSION
module load $PYTHON_MODULE

## Install Python libraries

# pip
#easy_install pip
PIP_PATH=$(which pip)
$PIP_PATH install --upgrade pip
 
# cython
$PIP_PATH install --upgrade cython

# docopt
$PIP_PATH install --upgrade docopt

# numpy
$PIP_PATH install --upgrade numpy

# biopython
$PIP_PATH install --upgrade biopython

# Markdown
$PIP_PATH install --upgrade markdown
$PIP_PATH install --upgrade pymdown-extensions

#Pands
$PIP_PATH install --upgrade pandas
$PIP_PATH install --upgrade tabulate

# python-dateutil
$PIP_PATH install --upgrade python-dateutil

# pyparsing
$PIP_PATH install --upgrade pyparsing

# matplotlib
$PIP_PATH install --upgrade matplotlib

# HTseq
$PIP_PATH install --upgrade HTseq

# bedtools-python has no version (version from doc: 0.1.0): install from master
$PIP_PATH install --upgrade  https://github.com/arq5x/bedtools-python/archive/master.zip

# PyVCF
$PIP_PATH install --upgrade PyVCF3

# pysam
$PIP_PATH install --upgrade pysam

# networkx
$PIP_PATH install --upgrade networkx

# sklearn
$PIP_PATH install --upgrade sklearn

#A permissions
chmod -R ug+rwX,o+rX-w $PYTHONHOME
