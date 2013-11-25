#!/bin/sh

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# You should probably also delete this commented-out header and the ## comments
################################################################################


#
# Software_name  ## TO BE MODIFIED WITH e.g. BLAST, HMMER, SAMtools, etc.
SOFTWARE="varscan2"
VERSION="2.3.6"  
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
wget http://downloads.sourceforge.net/project/varscan/VarScan.v$VERSION.jar  ## TO BE MODIFIED WITH SPECIFIC URL

# Add permissions and install software
chmod -R ug+rwX ./*
mv VarScan.v$VERSION.jar $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;  
}
module-whatis \"$SOFTWARE  \" ;  
                      
set             root                \$::env(MUGQIC_INSTALL_HOME_DEV)/software/$SOFTWARE ;
setenv          ${SOFTWARE}_JAR     \$root/VarScan.v$VERSION.jar ; 
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE
chmod -R ug+rwX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD


