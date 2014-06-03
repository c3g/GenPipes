#!/bin/bash

###################
################### GATK
###################
echo "You must download and install GATK manually";
echo "It can be found here:"
echo "http://www.broadinstitute.org/gatk/download"
VERSION=3.1-1
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/GenomeAnalysisTK/GenomeAnalysisTK-$VERSION
mkdir -p $INSTALL_PATH
echo "Install it here: $INSTALL_PATH"
echo "And keep the .tar.bz2 archive here: $MUGQIC_INSTALL_HOME/software/GenomeAnalysisTK/archive/"
chmod -R g+w $INSTALL_PATH

echo "#%Module1.0

proc ModulesHelp { } {
        puts stderr "\tadd  GenomeAnalysisTK"
}

module-whatis "The Broads toolsuite to work with resequencing"

set             root         $::env(MUGQIC_INSTALL_HOME)/software/GenomeAnalysisTK/GenomeAnalysisTK-$VERSION
setenv          GATK_JAR     \$root/GenomeAnalysisTK.jar
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/GenomeAnalysisTK
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/GenomeAnalysisTK/

echo "Module is installed here: $MUGQIC_INSTALL_HOME/modulefiles/mugqic/GenomeAnalysisTK"
