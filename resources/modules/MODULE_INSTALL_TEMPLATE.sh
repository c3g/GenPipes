#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

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

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV  ## TO BE MODIFIED IF NECESSARY

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE

# Create install directory with permissions if necessary
if [[ ! -d $INSTALL_DIR ]]
then
  mkdir $INSTALL_DIR
  chmod ug+rwX,o+rX-w $INSTALL_DIR
fi

INSTALL_DOWNLOAD=$INSTALL_DIR/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
ARCHIVE=$SOFTWARE-$VERSION.(zip|tar.gz|tar.bz2)  ## TO BE MODIFIED WITH SPECIFIC ARCHIVE
# If archive was previously downloaded, use the local one, otherwise get it from remote site
if [[ -f ${!INSTALL_HOME}/archive/$ARCHIVE ]]
then
  echo "Archive $ARCHIVE already in ${!INSTALL_HOME}/archive/: using it..."
  cp -a ${!INSTALL_HOME}/archive/$ARCHIVE .
else
  echo "Archive $ARCHIVE not in ${!INSTALL_HOME}/archive/: downloading it..."
  wget http://www.software_lab.org/download/$ARCHIVE -O $ARCHIVE  ## TO BE MODIFIED WITH SPECIFIC URL
fi
(unzip|tar zxvf|tar jxvf) $ARCHIVE  ## TO BE MODIFIED WITH SPECIFIC COMMAND

SOFTWARE_DIR=$SOFTWARE-$VERSION  ## TO BE MODIFIED WITH SPECIFIC SOFTWARE DIRECTORY IF NECESSARY
cd $SOFTWARE_DIR
./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR  ## TO BE ADDED AND MODIFIED IF NECESSARY
make  ## TO BE ADDED AND MODIFIED IF NECESSARY

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX,o+rX-w .
mv -i $SOFTWARE_DIR $INSTALL_DIR/
# Store archive if not already present or if different from the previous one
if [[ ! -f ${!INSTALL_HOME}/archive/$ARCHIVE || `diff ${!INSTALL_HOME}/archive/$ARCHIVE $ARCHIVE` ]]
then
  mv -i $ARCHIVE ${!INSTALL_HOME}/archive/
fi

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \" ;  ## TO BE MODIFIED WITH DETAILED HELP IF ANY
}
module-whatis \"$SOFTWARE\" ;  ## TO BE MODIFIED WITH DETAILED DESCRIPTION IF ANY

set             root                \$::env($INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin ;  ## TO BE ADDED IF NECESSARY
prepend-path    PATH                \$root/other_tools/bin ;  ## TO BE ADDED AND MODIFIED IF NECESSARY
setenv          ${SOFTWARE}_JAR     \$root/$SOFTWARE-$VERSION.jar ;  ## TO BE ADDED AND MODIFIED IF NECESSARY
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Set module directory path by lowercasing $INSTALL_HOME and removing '_install_home' in it
MODULE_DIR=${!INSTALL_HOME}/modulefiles/`echo ${INSTALL_HOME,,} | sed 's/_install_home//'`/$SOFTWARE

# Create module directory with permissions if necessary
if [[ ! -d $MODULE_DIR ]]
then
  mkdir $MODULE_DIR
  chmod ug+rwX,o+rX-w $MODULE_DIR
fi

# Add permissions and install module
chmod ug+rwX,o+rX-w $VERSION .version
mv $VERSION .version $MODULE_DIR/

# Clean up temporary installation files if any
cd
rm -rf $INSTALL_DOWNLOAD
