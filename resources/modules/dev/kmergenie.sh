#!/bin/bash

#
# KmerGenie
#

SOFTWARE="kmergenie"  ## TO BE MODIFIED WITH e.g. blast, hmmer, samtools, etc.
VERSION="1.6950"  ## TO BE MODIFIED WITH e.g. 2.2.28+, 3.0, 0.1.19, etc.
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
wget "http://kmergenie.bx.psu.edu/kmergenie-$VERSION.tar.gz" ## TO BE MODIFIED WITH SPECIFIC URL
tar zxvf $SOFTWARE-$VERSION.tar.gz  ## TO BE MODIFIED WITH SPECIFIC ARCHIVE
cd $SOFTWARE-$VERSION  ## TO BE MODIFIED WITH SPECIFIC $SOFTWARE-$VERSION IF DIFFERENT
source /etc/profile.d/modules.sh
module load mugqic_dev/R_Bioconductor/3.1.2_3.0
make  ## TO BE ADDED AND MODIFIED IF NECESSARY

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX .
mv -i $SOFTWARE-$VERSION $INSTALL_PATH  ## TO BE MODIFIED WITH SPECIFIC $SOFTWARE-$VERSION IF DIFFERENT
mv -i $SOFTWARE-$VERSION.tar.gz $MUGQIC_INSTALL_HOME_DEV/archive  ## TO BE MODIFIED WITH SPECIFIC $SOFTWARE-$VERSION.(zip|tar.gz|tar.bz2) IF DIFFERENT

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" 
}
module-whatis \"$SOFTWARE  \"  ## TO BE MODIFIED WITH DETAILED DESCRIPTION IF ANY

set             root                \$::env(MUGQIC_INSTALL_HOME_DEV)/software/$SOFTWARE/$SOFTWARE-$VERSION 
prepend-path    PATH                \$root
prepend-path    PATH                \$root/script  
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
