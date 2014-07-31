#!/bin/bash

# NOTE:
# This script is not fully portable because of BLAS and LAPACK. Pay attention to corresponding setenv paths below.
# TODO: we need to a way to keep track of the list of python packages + their versions for reporting, something a la help('modules')

###################
################### Python 
###################
VERSION=2.7.8
LIBVERSION=2.7
SUT_VERSION=5.4.1
# CYTHON_VERSION=0.20
# NUMPY_VERSION=1.8.0
# BIO_VERSION=1.63
# MATPLOT_VERSION=1.2.1
# HTSEQ_VERSION=0.5.4p1
PYVCF_VERSION=0.6.7

###MUGQIC_INSTALL_HOME=$SHARE_NOBACKUP/mugqic_prod => temporary mamouth transition
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/python/Python-${VERSION}
ARCHIVE_PATH=$MUGQIC_INSTALL_HOME/archive/python
PYTHONUSERBASE=${MUGQIC_INSTALL_HOME}/software/python/Python-${VERSION}
rm -rf $MUGQIC_INSTALL_HOME/software/python/Python-${VERSION}
mkdir -p $ARCHIVE_PATH $INSTALL_PATH

cd $MUGQIC_INSTALL_HOME/software/python/
wget http://www.python.org/ftp/python/$VERSION/Python-$VERSION.tgz
tar -xvf Python-$VERSION.tgz
cd Python-$VERSION
./configure --prefix=$INSTALL_PATH 
make -j8  
make install

cd ..
mv Python-$VERSION.tgz $ARCHIVE_PATH

# Module file 
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Adds Python $VERSION to your environment \"
}
module-whatis \"Adds Python $VERSION to your environment  \"

setenv  BLAS  			/software/libraries/GotoBLAS_LAPACK/shared/libblas.so
setenv  LAPACK		   /software/libraries/GotoBLAS_LAPACK/shared/liblapack.so

set             root               \$::env(MUGQIC_INSTALL_HOME)/software/python/Python-$VERSION
prepend-path    MANPATH            \$root/share/man              
prepend-path    PATH               \$root/bin
prepend-path    LD_LIBRARY_PATH    /software/libraries/GotoBLAS_LAPACK/shared:\$root/lib/
prepend-path    LIBRARY_PATH       /software/libraries/GotoBLAS_LAPACK/shared:\$root/lib/
prepend-path    CPATH              \$root/include:\$root/include/python${LIBVERSION}
prepend-path    PYTHONPATH      \$root/lib/python${LIBVERSION}/site-packages:\$root/lib/python${LIBVERSION}
#prepend-path    PYTHONHOME      \$root/lib
#setenv          PYTHONUSERBASE   \$root


" > $VERSION

#module-whatis	 Adds Python 2.6 to your environment  
#prepend-path	 MANPATH /software/tools/python-2.6.7/share/man 
#prepend-path	 PATH /software/tools/python-2.6.7/bin:/software/tools/swig-2.0.4/bin:/software/tools/wx-2.8.12/bin 
#prepend-path	 LD_LIBRARY_PATH /software/libraries/GotoBLAS_LAPACK/shared:/software/tools/python-2.6.7/lib:/software/tools/wx-2.8.12/lib 
#prepend-path	 CPATH /software/tools/python-2.6.7/include:/software/tools/wx-2.8.12/include 

#setenv		PYTHONPATH	\$root/lib
#setenv		PYTHONHOME	\$root/lib
# setenv                INTEL_LICENSE_FILE      $root/licenses
# setenv                FC                      gfortran
# setenv                F77                     gfortran
# setenv                CC                      gcc
#prepend-path    PATH               \$root/bin:/software/tools/swig-2.0.4/bin:/software/tools/wx-2.8.12/bin
#prepend-path    CPATH              $root/include:/software/tools/wx-2.8.12/include
#prepend-path    LD_LIBRARY_PATH    /software/libraries/GotoBLAS_LAPACK/shared:$root/lib:/software/tools/wx-2.8.12/lib




# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/python
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/python/

export LD_LIBRARY_PATH=/software/libraries/GotoBLAS_LAPACK/shared:${MUGQIC_INSTALL_HOME}/software/python/Python-${VERSION}/lib/python2.7:${LD_LIBRARY_PATH}
export LIBRARY_PATH=/software/libraries/GotoBLAS_LAPACK/shared:${MUGQIC_INSTALL_HOME}/software/python/Python-${VERSION}/lib/python2.7:${LIBRARY_PATH}
export CPATH=${MUGQIC_INSTALL_HOME}/software/python/Python-${VERSION}/include:${MUGQIC_INSTALL_HOME}/software/python/Python-${VERSION}/include/python${LIBVERSION}:${CPATH}
export PATH=${MUGQIC_INSTALL_HOME}/software/python/Python-${VERSION}/bin:${PATH}
export PYTHONPATH=${MUGQIC_INSTALL_HOME}/software/python/Python-${VERSION}/lib/python${LIBVERSION}/site-packages/:${PYTHONPATH}


#install setuptools => easy_install
mkdir -p $INSTALL_PATH/lib/python${LIBVERSION}/site-packages/
cd $INSTALL_PATH/lib/python${LIBVERSION}/site-packages/

wget --no-check-certificate https://pypi.python.org/packages/source/s/setuptools/setuptools-${SUT_VERSION}.tar.gz  
tar -xvf setuptools-${SUT_VERSION}.tar.gz 
cd setuptools-${SUT_VERSION}
python setup.py build
python setup.py install
cd ..
# 
# mv setuptools-${SUT_VERSION}.tar.gz $ARCHIVE_PATH

###################
################### Python Packages
###################
# Typically in e.g. /software/areas/genomics/software/python/Python-2.7.3/lib/python2.7/site-packages/
# mkdir -p $INSTALL_PATH/lib/python${LIBVERSION}/site-packages/
# cd $INSTALL_PATH/lib/python${LIBVERSION}/site-packages/
## Install cython (c code compiler for python)
easy_install Cython
# wget --no-check-certificate http://cython.org/release/Cython-${CYTHON_VERSION}.tar.gz
# tar -xvf Cython-${CYTHON_VERSION}.tar.gz
# cd Cython-${CYTHON_VERSION}
# python setup.py build  --prefix=$INSTALL_PATH
# python setup.py install --prefix=$INSTALL_PATH
# cd ..

## numpy
# module load mugqic/python/python-$VERSION
easy_install  numpy
# wget --no-check-certificate http://downloads.sourceforge.net/project/numpy/NumPy/${NUMPY_VERSION}/numpy-${NUMPY_VERSION}.tar.gz
# tar -xvf numpy-${NUMPY_VERSION}.tar.gz
# cd numpy-${NUMPY_VERSION}
# python setup.py build
# python setup.py install
# cd ..

##biopython
easy_install biopython
# wget --no-check-certificate http://biopython.org/DIST/biopython-${BIO_VERSION}.tar.gz
# tar -xvf biopython-${BIO_VERSION}.tar.gz
# cd biopython-${BIO_VERSION}
# python setup.py build
# python setup.py install
# cd ..

# ##toto
# ## matplotlib
# easy_install dateutil
easy_install matplotlib ### NOTE: didn't work on mamouth
# wget --no-check-certificate http://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-${MATPLOT_VERSION}/matplotlib-${MATPLOT_VERSION}.tar.gz
# tar -xvf matplotlib-${MATPLOT_VERSION}.tar.gz
# cd matplotlib-${MATPLOT_VERSION}
# python setup.py build # NOTE: no ssh -X from MacOS when building this
# python setup.py install
# cd ..

## Install HTSeq (HTseq is special, it creates an executalbe in python/bin)
easy_install HTSeq
# wget --no-check-certificate http://pypi.python.org/packages/source/H/HTSeq/HTSeq-${HTSEQ_VERSION}.tar.gz
# tar -xvf  HTSeq-${HTSEQ_VERSION}.tar.gz
# cd HTSeq-${HTSEQ_VERSION}
# python setup.py build
# python setup.py install
# cd ..

## Install bedtools-python
#no versionned repo (0.10)
easy_install https://github.com/arq5x/bedtools-python/archive/master.zip
# wget --no-check-certificate https://github.com/arq5x/bedtools-python/archive/master.zip
# gunzip master
# mv bedtools-python-master bedtools-python
# cd bedtools-python
# python setup.py build
# python setup.py install
# cd ..


## Install vcf
easy_install https://pypi.python.org/packages/source/P/PyVCF/PyVCF-${PYVCF_VERSION}.tar.gz
# wget --no-check-certificate https://pypi.python.org/packages/source/P/PyVCF/PyVCF-${PYVCF_VERSION}.tar.gz
# tar -xvf  PyVCF-${PYVCF_VERSION}.tar.gz
# cd PyVCF-${PYVCF_VERSION}
# python setup.py build
# python setup.py install
# cd ..

## matplotlib re do
easy_install http://labix.org/download/python-dateutil/python-dateutil-1.5.tar.gz
easy_install pyparsing
easy_install matplotlib

## RSeQC (easyinstall won't work)
# module load mugqic/python/2.7.5
# VERSION_RSeQC=2.3.9  
# wget https://downloads.sourceforge.net/project/rseqc/RSeQC-$VERSION_RSeQC.tar.gz
# tar xvf RSeQC-$VERSION_RSeQC.tar.gz
# cd RSeQC-$VERSION_RSeQC
# #module load mugqic/python/2.7.5 # wouldn't work with 2.7.8s
# python setup.py install
# cd ..
# rm -rf RSeQC*
# python -c 'from qcmodule import SAM'



#module-whatis	 HTSeq: Analysing high-throughput sequencing data with Python 
#prepend-path	 PATH /sb/programs/analyste/software/Python-2.7.3/bin/  # hmmm already done by python!!
# hmmm HT-Seq is just a python module and no two versions can co-exist. So more or less non-sensical to have a module.
# help('modules') would make a lot more sense
# Module file 
# echo "#%Module1.0
# proc ModulesHelp { } {
#        puts stderr \"\tMUGQIC - HTSeq dummy module \"
# module    load  mugqic/python
# }
# module-whatis \" HTSeq dummy module   \"
# " > $VERSION
# # version file
# echo "#%Module1.0
# set ModulesVersion \"$VERSION\"
# " > .version
# mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/htseq
# mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/htseq/
# import HTSeq to test

cd ..
rm -rf tmp
chmod 664 $MUGQIC_INSTALL_HOME/modulefiles/mugqic/python/$VERSION
chmod -R 775 $INSTALL_PATH

#####
## test install
#####
## module load mugqic/python/$VERSION
## htseq-count -h
## python
## from numpy import *
## from Bio import SeqIO
## from matplotlib import *
## from HTSeq import *
## from bedtools import *
## from vcf import *

