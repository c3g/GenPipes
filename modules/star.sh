#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#
# STAR
#

SOFTWARE=star
VERSION=2.4.0e

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME

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
ARCHIVE=${SOFTWARE^^}_$VERSION.tar.gz
# If archive was previously downloaded, use the local one, otherwise get it from remote site
if [[ -f ${!INSTALL_HOME}/archive/$ARCHIVE ]]
then
  echo "Archive $ARCHIVE already in ${!INSTALL_HOME}/archive/: using it..."
  cp -a ${!INSTALL_HOME}/archive/$ARCHIVE .
else
  echo "Archive $ARCHIVE not in ${!INSTALL_HOME}/archive/: downloading it..."
  wget http://github.com/alexdobin/STAR/archive/$ARCHIVE -O $ARCHIVE
fi
tar zxvf $ARCHIVE

SOFTWARE_DIR=${SOFTWARE^^}_$VERSION
# Remove "STAR-" prefix from top directory name
mv ${SOFTWARE^^}-$SOFTWARE_DIR $SOFTWARE_DIR
cd $SOFTWARE_DIR/source
make -j8

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
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                \$::env($INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE_DIR
prepend-path    PATH                \$root/source
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
