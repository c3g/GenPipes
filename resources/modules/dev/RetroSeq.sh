#!/bin/sh

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# You should probably also delete this commented-out header and the ## comments
################################################################################


#
# Software_name  ## TO BE MODIFIED WITH e.g. BLAST, HMMER, SAMtools, etc.
#

SOFTWARE=retroseq  ## TO BE MODIFIED WITH e.g. blast, hmmer, samtools, etc.
VERSION=0.1  ## TO BE MODIFIED WITH e.g. 2.2.28+, 3.0, 0.1.19, etc.
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/$SOFTWARE/$VERSION
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
git clone https://github.com/tk2/RetroSeq.git

# Add permissions and install software
chmod -R ug+rwX .
chmod -R o+rX .

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;  
}
module-whatis \"$SOFTWARE  \" ;  
                      
set             root                \$::env(MUGQIC_INSTALL_HOME_DEV)/software/$SOFTWARE/$VERSION/ ;  
prepend-path    PATH                \$root/RetroSeq/bin ;  
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE

# Clean up temporary installation files if any
#rm -rf $INSTALL_DOWNLOAD

