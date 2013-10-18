#!/bin/sh

#
# Tabix
#
SOFTWARE=tabix
VERSION=0.2.6
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, compile
wget http://downloads.sourceforge.net/project/samtools/$SOFTWARE/$SOFTWARE-$VERSION.tar.bz2
tar jxvf $SOFTWARE-$VERSION.tar.bz2
cd $SOFTWARE-$VERSION
make

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX .
mv -i $SOFTWARE-$VERSION $INSTALL_PATH  ## TO BE MODIFIED
mv -i $SOFTWARE-$VERSION.tar.bz2 $MUGQIC_INSTALL_HOME/archive  ## TO BE MODIFIED

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE add  Generic indexer for TAB-delimited genome position files\"
}
module-whatis \"$SOFTWARE Generic indexer for TAB-delimited genome position files \"
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION
prepend-path    PATH                \$root
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Install module
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
