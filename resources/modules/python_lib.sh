#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

# Python module where to install libs
PYTHON_VERSION=2.7.11
PYTHON_MODULE=mugqic/python/${PYTHON_VERSION}
module load $PYTHON_MODULE

## Install Python libraries

#pip
easy_install pip

# cython
#easy_install http://cython.org/release/Cython-0.23.4.tar.gz
pip install cython

# numpy
easy_install numpy
python -c 'import numpy; print numpy.__version__, numpy.__file__'

# biopython
easy_install biopython
python -c 'import Bio; print Bio.__version__, Bio.__file__'

easy_install http://labix.org/download/python-dateutil/python-dateutil-1.5.tar.gz
easy_install https://sourceforge.net/projects/pyparsing/files/pyparsing/pyparsing-2.1.1/pyparsing-2.1.1.tar.gz/download
easy_install https://github.com/matplotlib/matplotlib/archive/v1.5.1.tar.gz
python -c 'import matplotlib; print matplotlib.__version__, matplotlib.__file__'

# HTseq
easy_install https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.6.1p1.tar.gz
python -c 'import HTSeq; print HTSeq.__version__, HTSeq.__file__'

# bedtools-python has no version (version from doc: 0.1.0): install from master
easy_install https://github.com/arq5x/bedtools-python/archive/master.zip
python -c 'import bedtools; print  bedtools.__file__'

# PyVCF
easy_install https://pypi.python.org/packages/source/P/PyVCF/PyVCF-0.6.8.tar.gz
python -c 'import vcf; print vcf.__file__'

# pysam
easy_install https://github.com/pysam-developers/pysam/archive/v0.9.0.tar.gz
python -c 'import pysam; print pysam.__version__, pysam.__file__'

# nextworkx
easy_install https://github.com/networkx/networkx/archive/networkx-1.11.tar.gz
python -c 'import networkx; print networkx.__version__, networkx.__file__'


PIP_PATH=$(which pip)
# scikit-bio
${PIP_PATH} install scikit-bio
#futures
${PIP_PATH} install futures
# misopy
${PIP_PATH} install misopy

# Add permissions
chmod -R ug+rwX,o+rX-w $PYTHON_HOME

