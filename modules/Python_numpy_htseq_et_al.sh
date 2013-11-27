# NOTE:
# This script is not fully portable because of BLAS and LAPACK. Pay attention to corresponding setenv paths below.
# TODO: we need to a way to keep track of the list of python packages + their versions for reporting, something a la help('modules')

###################
################### Python 
###################
VERSION="2.7.6"
LIBVERSION="2.7"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/python/Python-${VERSION}
ARCHIVE_PATH=$MUGQIC_INSTALL_HOME/archive/python
mkdir -p $ARCHIVE_PATH $INSTALL_PATH/tmp/untar
cd $INSTALL_PATH/tmp
wget http://www.python.org/ftp/python/$VERSION/Python-$VERSION.tgz
tar -xvf Python-$VERSION.tgz -C untar
cd untar/Python-$VERSION
./configure --prefix=$INSTALL_PATH
make -j8
make install

cd ../..
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
prepend-path    LD_LIBRARY_PATH    /software/libraries/GotoBLAS_LAPACK/shared:\$root/lib
prepend-path    CPATH              \$root/include
prepend-path    PYTHONPATH      \$root/lib/python${LIBVERSION}
prepend-path    PYTHONHOME      \$root/lib
setenv          PYTHONUSERBASE   \$root


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

export PATH=${MUGQIC_INSTALL_HOME}/software/python/Python-${VERSION}/bin:${PATH}
export PYTHONPATH=${MUGQIC_INSTALL_HOME}/software/python/Python-${VERSION}/lib/python${LIBVERSION}:${PYTHONPATH}
export PYTHONUSERBASE=${MUGQIC_INSTALL_HOME}/software/python/Python-${VERSION}

#install setuptools => easy_install
SUT_VERSION=1.3.2
wget --no-check-certificate https://pypi.python.org/packages/source/s/setuptools/setuptools-${SUT_VERSION}.tar.gz  
tar -xvf setuptools-${SUT_VERSION}.tar.gz -C untar
cd untar/setuptools-${SUT_VERSION}
python setup.py build
python setup.py install
cd ../..

mv setuptools-${SUT_VERSION}.tar.gz $ARCHIVE_PATH

###################
################### Python Packages
###################
# Typically in e.g. /software/areas/genomics/software/python/Python-2.7.3/lib/python2.7/site-packages/
# mkdir -p $INSTALL_PATH/lib/python2.7/site-packages/
# cd $INSTALL_PATH/lib/python2.7/site-packages/

## Install cython (c code compiler for python)
easy_install cython

## numpy
easy_install numpy
# module load mugqic/python/python-$VERSION
# VERSION="1.8.0"
# wget http://downloads.sourceforge.net/project/numpy/NumPy/$VERSION/numpy-$VERSION.tar.gz
# tar -xvf numpy-$VERSION.tar.gz
# cd numpy-$VERSION
# python setup.py build
# python setup.py install
# cd ..

##biopython
easy_install biopython
# VERSION="1.62"
# wget http://biopython.org/DIST/biopython-${VERSION}.tar.gz
# tar -xvf biopython-${VERSION}.tar.gz
# cd biopython-${VERSION}
# python setup.py build
# python setup.py install
# cd ..

##toto
## matplotlib
easy_install matplotlib
# VERSION="1.2.1"
# wget "http://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-"$VERSION"/matplotlib-"$VERSION".tar.gz"
# tar -xvf matplotlib-$VERSION.tar.gz
# cd matplotlib-$VERSION
# python setup.py build # NOTE: no ssh -X from MacOS when building this
# python setup.py install
# cd ..

## Install HTSeq (HTseq is special, it creates an executalbe in python/bin)
easy_install HTSeq
# VERSION="0.5.4p1"
# wget http://pypi.python.org/packages/source/H/HTSeq/HTSeq-$VERSION.tar.gz
# tar -xvf  HTSeq-$VERSION.tar.gz
# cd HTSeq-$VERSION
# python setup.py build
# python setup.py install

## Install bedtools-python
easy_install https://github.com/arq5x/bedtools-python/archive/master.zip

## Install vcf
easy_install PyVCF

chmod 775 ${PYTHONUSERBASE}



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



