#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

# Python module where to install libs
PYTHON_VERSION=2.7.10_qiime
PYTHON_MODULE=mugqic/python/${PYTHON_VERSION}
module load $PYTHON_MODULE

## Install Python libraries

#pip
easy_install pip

# cython
easy_install cython

# numpy
easy_install numpy
python -c 'import numpy; print numpy.__version__, numpy.__file__'

# biopython
easy_install biopython
python -c 'import Bio; print Bio.__version__, Bio.__file__'

# matplotlib requires dateutil and pyparsing dependencies
easy_install python-dateutil
easy_install pyparsing
easy_install matplotlib==1.4.3
python -c 'import matplotlib; print matplotlib.__version__, matplotlib.__file__'

# HTseq
easy_install HTSeq
python -c 'import HTSeq; print HTSeq.__version__, HTSeq.__file__'

# bedtools-python has no version (version from doc: 0.1.0): install from master
easy_install https://github.com/arq5x/bedtools-python/archive/master.zip
python -c 'import bedtools; print  bedtools.__file__'

# PyVCF
easy_install PyVCF
python -c 'import vcf; print vcf.__file__'

# pysam
easy_install pysam
python -c 'import pysam; print pysam.__version__, pysam.__file__'

# nextworkx
easy_install networkx
python -c 'import networkx; print networkx.__version__, networkx.__file__'

PIP_PATH=$(which pip)
# TEToolkit
${PIP_PATH} install TEToolkit
# qiime
${PIP_PATH} install qiime

# Add permissions
chmod -R ug+rwX,o+rX-w $PYTHON_HOME
