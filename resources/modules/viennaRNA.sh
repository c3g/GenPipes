#!/bin/bash

SOFTWARE="ViennaRNA"  ## TO BE MODIFIED WITH e.g. blast, hmmer, samtools, etc.
VERSION="1.8.3"  
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
wget "http://www.tbi.univie.ac.at/~ivo/RNA/ViennaRNA-$VERSION.tar.gz"
tar -xvf ViennaRNA-$VERSION.tar.gz
cd $SOFTWARE-$VERSION 

# MANUAL FIX: sprintf fix
#MANUAL: nano RNAforester/src/rnafuncs.cpp  , add "#include <stdio.h>"
# http://missingreadme.wordpress.com/2010/11/08/how-to-install-the-vienna-rna-package/
./configure --prefix=$INSTALL_PATH/$SOFTWARE-$VERSION --datadir=$INSTALL_PATH/$SOFTWARE-$VERSION
# FIX: perl lib custom loc
cd Perl
perl Makefile.PL PREFIX=$MUGQIC_INSTALL_HOME/software/perl5libs
cd ..
make
make install

# Add permissions and install software
cd $INSTALL_PATH
rm -rf tmp
chmod -R ug+rwX  $INSTALL_PATH
chmod -R o+rX  $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;  
}
module-whatis \"$SOFTWARE  \" ; 
               
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION ; 
prepend-path    PATH                \$root/bin ;  
prepend-path    PATH                \$root/ViennaRNA/bin ;  
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

