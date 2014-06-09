#!/bin/bash

###################
################### BWA
###################
# tpx patch can be found here:
# ftp://ftp.conveysupport.com/outgoing/bwa/bwa-0.6.2-tpx.patch
VERSION="0.7.9a"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/bwa/
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download
wget http://downloads.sourceforge.net/project/bio-bwa/bwa-$VERSION.tar.bz2
tar xvjf bwa-$VERSION.tar.bz2

# Compile
cd bwa-$VERSION
#
# BWA
#

SOFTWARE=bwa
VERSION=0.7.9a

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME

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
ARCHIVE=$SOFTWARE-$VERSION.tar.bz2
# If archive was previously downloaded, use the local one, otherwise get it from remote site
if [[ -f ${!INSTALL_HOME}/archive/$ARCHIVE ]]
then
  echo "Archive $ARCHIVE already in ${!INSTALL_HOME}/archive/: using it..."
  cp -a ${!INSTALL_HOME}/archive/$ARCHIVE .
else
  echo "Archive $ARCHIVE not in ${!INSTALL_HOME}/archive/: downloading it..."
  wget http://downloads.sourceforge.net/project/bio-bwa/$ARCHIVE
fi
tar jxvf $ARCHIVE

SOFTWARE_DIR=$SOFTWARE-$VERSION
cd $SOFTWARE_DIR
>>>>>>> a441cda364f569f2c231fce63b9ed69270eaaec1
make -j8

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
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                \$::env($INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE_DIR
prepend-path    PATH                \$root
" > $VERSION

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
