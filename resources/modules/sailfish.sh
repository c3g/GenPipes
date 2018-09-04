#!/bin/bash

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# You should probably also delete this commented-out header and the ## comments
################################################################################

#
# Software_name  ## TO BE MODIFIED WITH e.g. BLAST, HMMER, SAMtools, etc.
#

SOFTWARE=sailfish  ## TO BE MODIFIED WITH e.g. blast, hmmer, samtools, etc.
VERSION=0.6.3  ## TO BE MODIFIED WITH e.g. 2.2.28+, 3.0, 0.1.19, etc.

# 'MUGQIC_INSTALL_HOME_TMP' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME_TMP  ## TO BE MODIFIED IF NECESSARY

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE

# Create install directory with permissions if necessary
if [[ ! -d $INSTALL_DIR ]]
then
  mkdir $INSTALL_DIR
  chmod ug+rwX,o+rX $INSTALL_DIR
fi

INSTALL_DOWNLOAD=$INSTALL_DIR/tmp
mkdir $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
ARCHIVE=Sailfish-$VERSION-Linux_x86-64.tar.gz  ## TO BE MODIFIED WITH SPECIFIC ARCHIVE

# 
# If archive was previously downloaded, use the local one, otherwise get it from remote site
if [[ -f ${!INSTALL_HOME}/archive/$ARCHIVE ]]
then
  echo "Archive $ARCHIVE already in ${!INSTALL_HOME}/archive/: using it..."
  cp -a ${!INSTALL_HOME}/archive/$ARCHIVE .
else
  echo "Archive $ARCHIVE not in ${!INSTALL_HOME}/archive/: downloading it..."
   wget  https://github.com/kingsfordgroup/sailfish/releases/download/v$VERSION/$ARCHIVE  -O $ARCHIVE
	
fi
tar xvf $ARCHIVE  



# https://github.com/kingsfordgroup/sailfish/blob/master/README.md
SOFTWARE_DIR=Sailfish-$VERSION-Linux_x86-64
# cd $SOFTWARE_DIR
# mkdir build
# cd build
#
# module load cmake/2.8.8 gcc/4.8.2  # mammouth
# #module load cmake/2.8.12.2 gcc/4.9.1  # , guillimin
# # abacus
# cmake -DFETCH_BOOST=TRUE ..
#
#
# make -j12
#
#
# make install



# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX,o+rX .
mv -i $SOFTWARE_DIR $INSTALL_DIR
# Store archive if not already present or if different from the previous one
if [[ ! -f ${!INSTALL_HOME}/archive/$ARCHIVE || `diff ${!INSTALL_HOME}/archive/$ARCHIVE $ARCHIVE` ]]
then
  mv -i $ARCHIVE ${!INSTALL_HOME}/archive/
fi

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \" ; 
}
module-whatis \"$SOFTWARE\" http://www.cs.cmu.edu/~ckingsf/software/sailfish/index.html ;  

set             root                \$::env($INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin ; 
prepend-path    LD_LIBRARY_PATH      \$root/lib ;  
" > $VERSION
# export =$PWD/lib/:$LD_LIBRARY_PATH

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Set module directory path by removing '_INSTALL_HOME' in $INSTALL_HOME and lowercasing the result
MODULE_DIR=${!INSTALL_HOME}/modulefiles/`echo ${INSTALL_HOME/_INSTALL_HOME/} | tr '[:upper:]' '[:lower:]'`/$SOFTWARE

# Create module directory with permissions if necessary
if [[ ! -d $MODULE_DIR ]]
then
  mkdir $MODULE_DIR
  chmod ug+rwX,o+rX $MODULE_DIR
fi

# Add permissions and install module
chmod ug+rwX,o+rX $VERSION .version
mv $VERSION .version $MODULE_DIR

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
