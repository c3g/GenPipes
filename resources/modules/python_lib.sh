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
easy_install http://cython.org/release/Cython-0.23.4.tar.gz

# numpy
easy_install http://downloads.sourceforge.net/project/numpy/NumPy/1.10.1/numpy-1.10.1.tar.gz
python -c 'import numpy; print numpy.__version__, numpy.__file__'

# biopython
easy_install http://biopython.org/DIST/biopython-1.66.tar.gz
python -c 'import Bio; print Bio.__version__, Bio.__file__'

# matplotlib requires dateutil and pyparsing dependencies (matplotlib is stuck at version 1.4.3 for qiime compatibility)
easy_install http://labix.org/download/python-dateutil/python-dateutil-1.5.tar.gz
easy_install http://downloads.sourceforge.net/project/pyparsing/pyparsing/pyparsing-2.0.6/pyparsing-2.0.6.tar.gz
easy_install https://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.4.3/matplotlib-1.4.3.tar.gz
python -c 'import matplotlib; print matplotlib.__version__, matplotlib.__file__'

# HTseq
easy_install https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.6.1p1.tar.gz
python -c 'import HTSeq; print HTSeq.__version__, HTSeq.__file__'

# bedtools-python has no version (version from doc: 0.1.0): install from master
easy_install https://github.com/arq5x/bedtools-python/archive/master.zip
python -c 'import bedtools; print  bedtools.__file__'

# PyVCF
easy_install https://pypi.python.org/packages/source/P/PyVCF/PyVCF-0.6.7.tar.gz
python -c 'import vcf; print vcf.__file__'

# pysam
easy_install https://github.com/pysam-developers/pysam/archive/v0.8.4.tar.gz
python -c 'import pysam; print pysam.__version__, pysam.__file__'

# nextworkx
easy_install http://networkx.lanl.gov/download/networkx/networkx-1.7.tar.gz
python -c 'import networkx; print networkx.__version__, networkx.__file__'


PIP_PATH=$(which pip)
# scikit-bio
${PIP_PATH} install scikit-bio
# qiime
${PIP_PATH} install qiime



# Add permissions
chmod -R ug+rwX,o+rX-w $PYTHON_HOME
