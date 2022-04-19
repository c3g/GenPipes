#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

# Python module where to install libs
PYTHON_VERSION=3.10.4
PYTHON_MODULE=mugqic/python/$PYTHON_VERSION
module load $PYTHON_MODULE

## Install Python libraries

# pip
#easy_install pip
PIP_PATH=$(which pip)
 
# cython
#easy_install http://cython.org/release/Cython-0.23.4.tar.gz
$PIP_PATH install --upgrade cython

# docopt
$PIP_PATH install --upgrade docopt

# numpy
$PIP_PATH install --upgrade numpy

# biopython
$PIP_PATH install --upgrade biopython

# python-dateutil
#easy_install http://labix.org/download/python-dateutil/python-dateutil-1.5.tar.gz
$PIP_PATH install --upgrade python-dateutil

# pyparsing
#easy_install https://sourceforge.net/projects/pyparsing/files/pyparsing/pyparsing-2.1.1/pyparsing-2.1.1.tar.gz/download
$PIP_PATH install --upgrade pyparsing

# matplotlib
#easy_install https://github.com/matplotlib/matplotlib/archive/v1.5.1.tar.gz
$PIP_PATH install --upgrade matplotlib

# HTseq
#easy_install https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.6.1p1.tar.gz
$PIP_PATH install --upgrade HTseq

# bedtools-python has no version (version from doc: 0.1.0): install from master
$PIP_PATH install --upgrade  https://github.com/arq5x/bedtools-python/archive/master.zip

# PyVCF
#easy_install https://pypi.python.org/packages/source/P/PyVCF/PyVCF-0.6.8.tar.gz
$PIP_PATH install --upgrade PyVCF3

# pysam
#easy_install https://github.com/pysam-developers/pysam/archive/v0.9.0.tar.gz
$PIP_PATH install --upgrade pysam

# networkx
$PIP_PATH install --upgrade networkx

# futures
$PIP_PATH install --upgrade futures
# future
$PIP_PATH install --upgrade future

# misopy
$PIP_PATH install --upgrade misopy

# qiime
$PIP_PATH install --upgrade qiime

# TEToolKit
$PIP_PATH install --upgrade TEToolkit

# deepTools
$PIP_PATH install --upgrade deeptools

# sklearn
$PIP_PATH install --upgrade sklearn

# RSeQC
$PIP_PATH install --upgrade RSeQC

#A permissions
chmod -R ug+rwX,o+rX-w $PYTHONHOME
