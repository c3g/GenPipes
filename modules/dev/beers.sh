#!/bin/bash

#
# Software_name  ## TO BE MODIFIED WITH e.g. BLAST, HMMER, SAMtools, etc.
#
SOFTWARE=beers  ## TO BE MODIFIED WITH e.g. blast, hmmer, samtools, etc.
VERSION=0.0.0  ## TO BE MODIFIED WITH e.g. 2.2.28+, 3.0, 0.1.19, etc.
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE/$SOFTWARE-$VERSION
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
wget -q http://www.cbil.upenn.edu/BEERS/beers.tar 
tar -xvf beers.tar
mv  beers.tar $MUGQIC_INSTALL_HOME/archive

# Get the template config files
mkdir -p config_refseq && cd config_refseq
wget -q http://itmat.rum.s3.amazonaws.com/simulator_config_refseq.tar.gz
tar -xvf simulator_*
rm simulator_*.tar.gz
cd ..

# Install necessary perl module
module load mugqic/perl
cpan Math::Random

# Add permissions and install software
chmod -R ug+rwX .
chmod -R o+rX .

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE RNA-seq reads simulator \" ;  
}
module-whatis \"$SOFTWARE  \" ;  
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION ;  
setenv          BEERS_HOME     \$root ; 
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE


