#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

# Python module where to install libs
PYTHON_VERSION=2.7.13
PYTHON_MODULE=mugqic/python/$PYTHON_VERSION
module load $PYTHON_MODULE

## Install Python libraries

# pip
easy_install pip
PIP_PATH=$(which pip)
 
# cython
#easy_install http://cython.org/release/Cython-0.23.4.tar.gz
$PIP_PATH install --upgrade cython
python -c 'import cython; print cython.__version__, cython.__file__'

# docopt
$PIP_PATH install --upgrade docopt
python -c 'import docopt; print docopt.__version__, docopt.__file__'

# numpy
$PIP_PATH install --upgrade numpy
python -c 'import numpy; print numpy.__version__, numpy.__file__'

# biopython
$PIP_PATH install --upgrade biopython
python -c 'import Bio; print Bio.__version__, Bio.__file__'

# python-dateutil
#easy_install http://labix.org/download/python-dateutil/python-dateutil-1.5.tar.gz
$PIP_PATH install --upgrade python-dateutil
python -c 'import dateutil; print dateutil.__version__, dateutil.__file__'

# pyparsing
#easy_install https://sourceforge.net/projects/pyparsing/files/pyparsing/pyparsing-2.1.1/pyparsing-2.1.1.tar.gz/download
$PIP_PATH install --upgrade pyparsing
python -c 'import pyparsing; print pyparsing.__version__, pyparsing.__file__'

# matplotlib
#easy_install https://github.com/matplotlib/matplotlib/archive/v1.5.1.tar.gz
$PIP_PATH install --upgrade matplotlib
python -c 'import matplotlib; print matplotlib.__version__, matplotlib.__file__'

# HTseq
#easy_install https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.6.1p1.tar.gz
$PIP_PATH install --upgrade HTseq
python -c 'import HTSeq; print HTSeq.__version__, HTSeq.__file__'

# bedtools-python has no version (version from doc: 0.1.0): install from master
$PIP_PATH install --upgrade  https://github.com/arq5x/bedtools-python/archive/master.zip
python -c 'import bedtools; print bedtools.__file__'

# PyVCF
#easy_install https://pypi.python.org/packages/source/P/PyVCF/PyVCF-0.6.8.tar.gz
$PIP_PATH install --upgrade PyVCF
python -c 'import vcf; print vcf.__file__'

# pysam
#easy_install https://github.com/pysam-developers/pysam/archive/v0.9.0.tar.gz
$PIP_PATH install --upgrade pysam
python -c 'import pysam; print pysam.__version__, pysam.__file__'

# networkx
$PIP_PATH install --upgrade networkx
python -c 'import networkx; print networkx.__version__, networkx.__file__'

# futures
$PIP_PATH install --upgrade futures
# future
$PIP_PATH install --upgrade future
python -c 'import future; print future.__version__, future.__file__'

# misopy
$PIP_PATH install --upgrade misopy
python -c 'import misopy; print misopy.__version__, misopy.__file__'

# qiime
$PIP_PATH install --upgrade qiime
python -c 'import qiime; print qiime.__version__, qiime.__file__'

# TEToolKit
$PIP_PATH install --upgrade TEToolkit
python -c 'import TEToolkit; print TEToolkit.__file__'

# deepTools
$PIP_PATH install --upgrade deeptools
python -c 'import deeptools; print deeptools.__file__'

# sklearn
$PIP_PATH install --upgrade sklearn
python -c 'import sklearn; print sklearn.__version__; print sklearn.__file__'

# RSeQC
$PIP_PATH install --upgrade RSeQC

#A permissions
chmod -R ug+rwX,o+rX-w $PYTHONHOME
