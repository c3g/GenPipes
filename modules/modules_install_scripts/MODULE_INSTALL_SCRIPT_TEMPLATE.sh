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

SOFTWARE=software_name  ## TO BE MODIFIED WITH e.g. blast, hmmer, samtools, etc.
VERSION=0.0.0  ## TO BE MODIFIED WITH e.g. 2.2.28+, 3.0, 0.1.19, etc.
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
wget http://www.software_lab.org/download/$SOFTWARE-$VERSION.(zip|tar.gz|tar.bz2)  ## TO BE MODIFIED WITH SPECIFIC URL
(unzip|tar zxvf|tar jxvf) $SOFTWARE-$VERSION.(zip|tar.gz|tar.bz2)  ## TO BE MODIFIED WITH SPECIFIC ARCHIVE
cd $SOFTWARE-$VERSION  ## TO BE MODIFIED WITH SPECIFIC $SOFTWARE-$VERSION IF DIFFERENT
./configure --prefix=$INSTALL_PATH/$SOFTWARE-$VERSION  ## TO BE ADDED AND MODIFIED IF NECESSARY
make  ## TO BE ADDED AND MODIFIED IF NECESSARY

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX .
mv -i $SOFTWARE-$VERSION $INSTALL_PATH  ## TO BE MODIFIED WITH SPECIFIC $SOFTWARE-$VERSION IF DIFFERENT
mv -i $SOFTWARE-$VERSION.(zip|tar.gz|tar.bz2) $MUGQIC_INSTALL_HOME/archive  ## TO BE MODIFIED WITH SPECIFIC $SOFTWARE-$VERSION.(zip|tar.gz|tar.bz2) IF DIFFERENT

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;  ## TO BE MODIFIED WITH DETAILED DESCRIPTION IF ANY
}
module-whatis \"$SOFTWARE  \" ;  ## TO BE MODIFIED WITH DETAILED DESCRIPTION IF ANY
                      
prereq          mugqic/python/2.7.3 ;  ## TO BE ADDED WITH SPECIFIC PREREQUISITE IF NECESSARY
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION ;  ## TO BE MODIFIED WITH SPECIFIC $SOFTWARE-$VERSION IF DIFFERENT
prepend-path    PATH                \$root/bin ;  ## TO BE ADDED IF NECESSARY
prepend-path    PATH                \$root/other_tools/bin ;  ## TO BE ADDED AND MODIFIED IF NECESSARY
setenv          ${SOFTWARE}_JAR     \$root/$SOFTWARE-$VERSION.jar ;  ## TO BE ADDED AND MODIFIED IF NECESSARY
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE
chmod -R ug+rwX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
