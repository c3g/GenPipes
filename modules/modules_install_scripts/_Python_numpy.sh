

###################
################### Python
###################
VERSION="2.7.3"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/python/Python-$VERSION
mkdir -p $INSTALL_PATH
wget "http://www.python.org/ftp/python/$VERSION/Python-$VERSION.tgz"
tar -xvf Python-$VERSION.tgz
cd Python-$VERSION
./configure --prefix=$INSTALL_PATH
make -j8
make install


# Module file (TODO: BLAS and LAPACK are not portable)
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
prepend-path    LD_LIBRARY_PATH    \$root/lib
prepend-path    CPATH              \$root/include

#setenv		PYTHONPATH	\$root/lib
#setenv		PYTHONHOME	\$root/lib
# setenv                INTEL_LICENSE_FILE      $root/licenses
# setenv                FC                      gfortran
# setenv                F77                     gfortran
# setenv                CC                      gcc
#prepend-path    PATH               \$root/bin:/software/tools/swig-2.0.4/bin:/software/tools/wx-2.8.12/bin
#prepend-path    CPATH              $root/include:/software/tools/wx-2.8.12/include
#prepend-path    LD_LIBRARY_PATH    /software/libraries/GotoBLAS_LAPACK/shared:$root/lib:/software/tools/wx-2.8.12/lib


" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/python
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/python/


## Install numpy
module load mugqic/python



## Install HTSeq
module load mugqic/python






