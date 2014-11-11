#!/bin/bash
###################
################### Beagle
###################
VERSION="4.r1128"
VERSION_DOC="4preview_27Jun13"

INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/beagle/beagle-${VERSION}
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
# Download
wget --content-disposition http://faculty.washington.edu/browning/beagle/b${VERSION}.jar
wget --content-disposition http://faculty.washington.edu/browning/beagle/b${VERSION_DOC}.pdf
chmod -R g+w $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Java tool to phase genomes \"
}
module-whatis \"MUGQIC - beagle \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/beagle/beagle-${VERSION}
setenv          BEAGLE_HOME        \$root
" > $VERSION
# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/beagle
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/beagle/
