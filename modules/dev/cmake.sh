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

SOFTWARE=cmake  ## TO BE MODIFIED WITH e.g. blast, hmmer, samtools, etc.
VERSION=2.8.12.1-Linux-i386  ## TO BE MODIFIED WITH e.g. 2.2.28+, 3.0, 0.1.19, etc.
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/$SOFTWARE/$VERSION
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
wget http://www.cmake.org/files/v2.8/$SOFTWARE-$VERSION.sh  
chmod 775 $SOFTWARE-$VERSION.sh

#you need answer : y and n
./$SOFTWARE-$VERSION.sh 

chmod -R ug+rwX .
chmod -R o+rX .
mv $SOFTWARE-$VERSION.sh $MUGQIC_INSTALL_HOME_DEV/archive  ## TO BE MODIFIED WITH SPECIFIC $SOFTWARE-$VERSION.(zip|tar.gz|tar.bz2) IF DIFFERENT

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;  
}
module-whatis \"$SOFTWARE : to build makefile \" ;
                      
set             root                \$::env(MUGQIC_INSTALL_HOME_DEV)/software/$SOFTWARE/$VERSION
prepend-path    PATH                \$root/bin ;
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
