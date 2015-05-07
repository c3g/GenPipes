#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

# Python module where to install libs
PYTHON_MODULE=mugqic/python/2.7.8
module load $PYTHON_MODULE

## Install Python libraries

# cython
easy_install http://cython.org/release/Cython-0.22.tar.gz

# numpy
easy_install http://downloads.sourceforge.net/project/numpy/NumPy/1.9.1/numpy-1.9.1.tar.gz
python -c 'import numpy; print numpy.__version__, numpy.__file__'

# biopython
easy_install http://biopython.org/DIST/biopython-1.65.tar.gz
python -c 'import Bio; print Bio.__version__, Bio.__file__'

# matplotlib requires dateutil and pyparsing dependencies
easy_install http://labix.org/download/python-dateutil/python-dateutil-1.5.tar.gz
easy_install http://downloads.sourceforge.net/project/pyparsing/pyparsing/pyparsing-2.0.3/pyparsing-2.0.3.tar.gz
easy_install https://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.4.2/matplotlib-1.4.2.tar.gz
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
easy_install https://github.com/pysam-developers/pysam/archive/v0.8.2.tar.gz
python -c 'import pysam; print pysam.__version__, pysam.__file__'

# nextworkx
easy_install http://networkx.lanl.gov/download/networkx/networkx-1.7.tar.gz
python -c 'import networkx; print networkx.__version__, networkx.__file__'

# Add permissions
chmod -R ug+rwX,o+rX-w $PYTHON_HOME
