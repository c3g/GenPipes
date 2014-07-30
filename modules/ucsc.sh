#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#
# UCSC 'kent' bioinformatics utilities
#

SOFTWARE=ucsc
# By default, the latest remote version will be downloaded and the version date set appropriately.
# To use a local archive specific version, uncomment and update VERSION
VERSION=20140212

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
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:

# Momentarily unset variable expansion error to test $VERSION
set +u
if [[ -z $VERSION ]]
then
  # If VERSION is not set, download the latest remote archive
  echo "Archive VERSION not set: downloading latest remote archive..."
  REMOTE_ARCHIVE=userApps.src.tgz
  wget http://hgdownload.cse.ucsc.edu/admin/exe/$REMOTE_ARCHIVE
  # Set VERSION with the archive last modification date
  VERSION=`stat --printf=%y $REMOTE_ARCHIVE | perl -pe 's/^(\d+)-(\d+)-(\d+).*/\1\2\3/'`
  ARCHIVE=$SOFTWARE-userApps-$VERSION.src.tgz
  mv $REMOTE_ARCHIVE $ARCHIVE
else
  # Try to find specific version in archive repository
  ARCHIVE=$SOFTWARE-userApps-$VERSION.src.tgz
  if [[ -f ${!INSTALL_HOME}/archive/$ARCHIVE ]]
  then
    echo "Archive $ARCHIVE already in ${!INSTALL_HOME}/archive/: using it..."
    cp -a ${!INSTALL_HOME}/archive/$ARCHIVE .
  else
    echo "Error: archive $ARCHIVE not found in ${!INSTALL_HOME}/archive/!"
    echo "Comment VERSION variable to download the latest remote version"
    exit 1
  fi
fi
set -u

tar zxvf $ARCHIVE

SOFTWARE_DIR=$SOFTWARE-$VERSION
mv userApps $SOFTWARE_DIR
cd $SOFTWARE_DIR
make

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX,o+rX .
mv -i $SOFTWARE_DIR $INSTALL_DIR/
# Store archive if not already present or if different from the previous one
if [[ ! -f ${!INSTALL_HOME}/archive/$ARCHIVE || `diff -q ${!INSTALL_HOME}/archive/$ARCHIVE $ARCHIVE` ]]
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
prepend-path    PATH                \$root/bin
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
mv $VERSION .version $MODULE_DIR/

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
