#!/bin/bash

###################
################### BWA
###################
# tpx patch can be found here:
# ftp://ftp.conveysupport.com/outgoing/bwa/bwa-0.6.2-tpx.patch
VERSION="20170322"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/parallel/
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download
wget http://ftp.gnu.org/gnu/parallel/parallel-${VERSION}.tar.bz2
tar xvjf parallel-$VERSION.tar.bz2
mv parallel-$VERSION parallel-${VERSION}-src
# Compile
cd parallel-${VERSION}-src
./configure --prefix=$MUGQIC_INSTALL_HOME/software/parallel/parallel-${VERSION}
make -j12
make install
cd ..
rm -r parallel-${VERSION}-src

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - GNU Parallel \"
}
module-whatis \"MUGQIC - GNU Parallel \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/parallel/parallel-${VERSION}
prepend-path    PATH               \$root/bin
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/parallel
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/parallel/

