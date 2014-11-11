#!/bin/bash

###################
################### MuTect
###################
echo "New versions of mutect must be built from tag releases";
echo "It can be found here:"
echo "https://github.com/broadinstitute/mutect/tree/master"
echo "Compile using the exact procedure (the same gatk version, and don't forget the git checkout tag <VERSION>";
VERSION=1.1.6
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/mutect/muTect-$VERSION

mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
unzip ${MUGQIC_INSTALL_HOME}/archive/muTect-1.1.6-bin.zip
chmod -R g+w $INSTALL_PATH

echo "#%Module1.0

proc ModulesHelp { } {
        puts stderr "\tadd MuTect"
}

module-whatis "The MuTect somatic caller"

set             root         $::env(MUGQIC_INSTALL_HOME)/software/mutect/muTect-$VERSION
setenv          MUTECT_JAR  \$root/muTect-1.1.6.jar
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/mutect
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/mutect/

echo "Module is installed here: $MUGQIC_INSTALL_HOME/modulefiles/mugqic/mutect"
