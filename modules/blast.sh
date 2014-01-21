#!/bin/sh

#
# BLAST
#

SOFTWARE=blast
VERSION=2.2.29+
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${VERSION%+}/ncbi-$SOFTWARE-$VERSION-x64-linux.tar.gz
tar zxvf ncbi-$SOFTWARE-$VERSION-x64-linux.tar.gz

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX .
chmod -R o+rX .
mv -i ncbi-$SOFTWARE-$VERSION $INSTALL_PATH
mv -i ncbi-$SOFTWARE-$VERSION-x64-linux.tar.gz $MUGQIC_INSTALL_HOME/archive

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - NCBI $SOFTWARE \" ;
}
module-whatis \"MUGQIC NCBI Blast: blast alignment \" ;
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/ncbi-$SOFTWARE-$VERSION ;
prepend-path    PATH                \$root/bin ;
prepend-path    BLASTDB             \$::env(MUGQIC_INSTALL_HOME)/genomes/blast_db ;
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

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
