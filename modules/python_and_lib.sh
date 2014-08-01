#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#
# Python
#

# NOTE:
# This script is not fully portable because of BLAS and LAPACK. Pay attention to corresponding setenv paths below.
# TODO: we need to a way to keep track of the list of python packages + their versions for reporting, something a la help('modules')

SOFTWARE=python
VERSION=2.7.8
# Remove the version last number
LIBVERSION=${VERSION%.[0-9]*}
SETUPTOOLS_VERSION=5.4.1
# CYTHON_VERSION=0.20
# NUMPY_VERSION=1.8.0
# BIOPYTHON_VERSION=1.63
# MATPLOTLIB_VERSION=1.2.1
# HTSEQ_VERSION=0.5.4p1
PYVCF_VERSION=0.6.7

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE

# Create install directory with permissions if necessary
if [[ ! -d $INSTALL_DIR ]]
then
  mkdir $INSTALL_DIR
  chmod ug+rwX,o+rX $INSTALL_DIR
fi

INSTALL_DOWNLOAD=$INSTALL_DIR/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Uppercase first P in python
ARCHIVE=${SOFTWARE^}-$VERSION.tgz
# If archive was previously downloaded, use the local one, otherwise get it from remote site
if [[ -f ${!INSTALL_HOME}/archive/$ARCHIVE ]]
then
  echo "Archive $ARCHIVE already in ${!INSTALL_HOME}/archive/: using it..."
  cp -a ${!INSTALL_HOME}/archive/$ARCHIVE .
else
  echo "Archive $ARCHIVE not in ${!INSTALL_HOME}/archive/: downloading it..."
  wget --no-check-certificate http://www.python.org/ftp/$SOFTWARE/$VERSION/$ARCHIVE
fi
tar zxvf $ARCHIVE

# Uppercase first P in python
SOFTWARE_DIR=${SOFTWARE^}-$VERSION
SOFTWARE_INSTALL_DIR=$INSTALL_DIR/$SOFTWARE_DIR
cd $SOFTWARE_DIR
./configure --prefix=$SOFTWARE_INSTALL_DIR
make -j8
rm -rf $SOFTWARE_INSTALL_DIR
make install

export LD_LIBRARY_PATH=/software/libraries/GotoBLAS_LAPACK/shared:$SOFTWARE_INSTALL_DIR/lib/python$LIBVERSION
export LIBRARY_PATH=/software/libraries/GotoBLAS_LAPACK/shared:$SOFTWARE_INSTALL_DIR/lib/python$LIBVERSION
export CPATH=$SOFTWARE_INSTALL_DIR/include:$SOFTWARE_INSTALL_DIR/include/python$LIBVERSION
export PATH=$SOFTWARE_INSTALL_DIR/bin:$PATH
export PYTHONPATH=$SOFTWARE_INSTALL_DIR/lib/python$LIBVERSION/site-packages

# Install setuptools => easy_install
mkdir -p $SOFTWARE_INSTALL_DIR/lib/python$LIBVERSION/site-packages
cd $SOFTWARE_INSTALL_DIR/lib/python$LIBVERSION/site-packages

wget --no-check-certificate https://pypi.python.org/packages/source/s/setuptools/setuptools-$SETUPTOOLS_VERSION.tar.gz
tar zxvf setuptools-$SETUPTOOLS_VERSION.tar.gz
cd setuptools-$SETUPTOOLS_VERSION
python setup.py build
python setup.py install


#
# Python packages
#

# Cython (c code compiler for python)
easy_install Cython
# cd $SOFTWARE_INSTALL_DIR/lib/python$LIBVERSION/site-packages
# wget --no-check-certificate http://cython.org/release/Cython-${CYTHON_VERSION}.tar.gz
# tar -xvf Cython-${CYTHON_VERSION}.tar.gz
# cd Cython-${CYTHON_VERSION}
# python setup.py build --prefix=$INSTALL_PATH
# python setup.py install --prefix=$INSTALL_PATH

# numpy
easy_install  numpy
# cd $SOFTWARE_INSTALL_DIR/lib/python$LIBVERSION/site-packages
# wget --no-check-certificate http://downloads.sourceforge.net/project/numpy/NumPy/${NUMPY_VERSION}/numpy-${NUMPY_VERSION}.tar.gz
# tar -xvf numpy-${NUMPY_VERSION}.tar.gz
# cd numpy-${NUMPY_VERSION}
# python setup.py build
# python setup.py install

# biopython
easy_install biopython
# cd $SOFTWARE_INSTALL_DIR/lib/python$LIBVERSION/site-packages
# wget --no-check-certificate http://biopython.org/DIST/biopython-${BIOPYTHON_VERSION}.tar.gz
# tar -xvf biopython-${BIOPYTHON_VERSION}.tar.gz
# cd biopython-${BIOPYTHON_VERSION}
# python setup.py build
# python setup.py install

# HTSeq (HTseq is special, it creates an executable in python/bin)
easy_install HTSeq
# cd $SOFTWARE_INSTALL_DIR/lib/python$LIBVERSION/site-packages
# wget --no-check-certificate http://pypi.python.org/packages/source/H/HTSeq/HTSeq-${HTSEQ_VERSION}.tar.gz
# tar -xvf  HTSeq-${HTSEQ_VERSION}.tar.gz
# cd HTSeq-${HTSEQ_VERSION}
# python setup.py build
# python setup.py install

# bedtools-python
# No versioned repo (0.10)
easy_install https://github.com/arq5x/bedtools-python/archive/master.zip
# cd $SOFTWARE_INSTALL_DIR/lib/python$LIBVERSION/site-packages
# wget --no-check-certificate https://github.com/arq5x/bedtools-python/archive/master.zip
# gunzip master
# mv bedtools-python-master bedtools-python
# cd bedtools-python
# python setup.py build
# python setup.py install

# PyVCF
easy_install https://pypi.python.org/packages/source/P/PyVCF/PyVCF-${PYVCF_VERSION}.tar.gz
# cd $SOFTWARE_INSTALL_DIR/lib/python$LIBVERSION/site-packages
# wget --no-check-certificate https://pypi.python.org/packages/source/P/PyVCF/PyVCF-${PYVCF_VERSION}.tar.gz
# tar -xvf  PyVCF-${PYVCF_VERSION}.tar.gz
# cd PyVCF-${PYVCF_VERSION}
# python setup.py build
# python setup.py install

# matplotlib
easy_install http://labix.org/download/python-dateutil/python-dateutil-1.5.tar.gz
easy_install pyparsing
# matplotlib first install fails: repeat it a second time to succeed (?!: pb with tornado 4.0 to be investigated)
set +e
easy_install matplotlib
easy_install matplotlib
set -e
# cd $SOFTWARE_INSTALL_DIR/lib/python$LIBVERSION/site-packages
# wget --no-check-certificate http://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-${MATPLOTLIB_VERSION}/matplotlib-${MATPLOTLIB_VERSION}.tar.gz
# tar -xvf matplotlib-${MATPLOTLIB_VERSION}.tar.gz
# cd matplotlib-${MATPLOTLIB_VERSION}
# python setup.py build # NOTE: no ssh -X from MacOS when building this
# python setup.py install



# Add permissions
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX,o+rX . $INSTALL_DIR/$SOFTWARE_DIR
# Store archive if not already present or if different from the previous one
if [[ ! -f ${!INSTALL_HOME}/archive/$ARCHIVE || `diff ${!INSTALL_HOME}/archive/$ARCHIVE $ARCHIVE` ]]
then
  mv -i $ARCHIVE ${!INSTALL_HOME}/archive/
fi

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

setenv  BLAS                        /software/libraries/GotoBLAS_LAPACK/shared/libblas.so
setenv  LAPACK                      /software/libraries/GotoBLAS_LAPACK/shared/liblapack.so

set             root                \$::env($INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
prepend-path    MANPATH             \$root/share/man
prepend-path    LIBRARY_PATH        /software/libraries/GotoBLAS_LAPACK/shared:\$root/lib/
prepend-path    LD_LIBRARY_PATH     /software/libraries/GotoBLAS_LAPACK/shared:\$root/lib/
prepend-path    CPATH               \$root/include:\$root/include/python${LIBVERSION}
prepend-path    PYTHONPATH          \$root/lib/python${LIBVERSION}/site-packages:\$root/lib/python${LIBVERSION}
" > $VERSION

################################################################################
# Everything between those lines should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Set module directory path by lowercasing $INSTALL_HOME and removing '_install_home' in it
MODULE_DIR=${!INSTALL_HOME}/modulefiles/`echo ${INSTALL_HOME,,} | sed 's/_install_home//'`/$SOFTWARE

# Create module directory with permissions if necessary
if [[ ! -d $MODULE_DIR ]]
then
  mkdir $MODULE_DIR
  chmod ug+rwX,o+rX $MODULE_DIR
fi

# Add permissions and install module
chmod ug+rwX,o+rX $VERSION .version
mv $VERSION .version $MODULE_DIR/

# Clean up temporary installation files if any
cd
rm -rf $INSTALL_DOWNLOAD

################################################################################

# Test module install
echo `seq -s - 80 | sed 's/[0-9]//g'`
echo "Testing $MODULE_DIR/$VERSION install..."
module load $MODULE_DIR/$VERSION
htseq-count -h
python -c '\
from numpy import *
from Bio import SeqIO
from matplotlib import *
from HTSeq import *
from bedtools import *
from vcf import *'
