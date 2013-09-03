#!/bin/bash
###################
################### beagle
###################
VERSION="3.3.2"
VERSION_DOC="${VERSION}_31Oct11"

INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/beagle/beagle-${VERSION}
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
# Download
wget --content-disposition http://faculty.washington.edu/browning/beagle/beagle.jar
wget --content-disposition http://faculty.washington.edu/browning/beagle/beagle_${VERSION_DOC}.pdf

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Java tool to phase genomes \"
}
module-whatis \"MUGQIC - beagle \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/beagle/beagle-${VERSION}
setenv          BEAGLE_HOME        \$root
" > $VERSION
# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/beagle
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/beagle/
