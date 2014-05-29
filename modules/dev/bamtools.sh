#!/bin/sh

###################
################### bamtools
###################
module load cmake
VERSION="2.3.0"
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/bamtools/bamtools-$VERSION
mkdir -p $INSTALL_PATH

git clone git://github.com/pezmaster31/bamtools.git
cd bamtools
mkdir -p build
cd build
cmake ..
make
cd ..
cp -r bin $INSTALL_PATH/

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Bamtools \"
}
module-whatis \"Ray Parallel genome assemblies for parallel DNA sequencing \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/bamtools/bamtools-$VERSION
prepend-path    PATH               \$root/bin
prepend-path    LD_LIBRARY_PATH    \$root/lib
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/bamtools
mv .version $VERSION $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/bamtools/


