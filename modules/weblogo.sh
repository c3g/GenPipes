#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#
# WebLogo
#

# WebLogo requires a recent version of ghostscript to create PNG and PDF output, and pdf2svg to generate SVG output.
# WebLogo version 3 is written in python. It is necessary to have python 2.5, 2.6 or 2.7 and the extension package numpy installed before WebLogo will run. 

SOFTWARE=weblogo
#VERSION=2.8.2
VERSION=3.3
# If WebLogo version >= 3, specify in which python version it must be intalled
PYTHON_VERSION=2.7.8

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
if [[ $VERSION == "2.8.2" ]]
then
  ARCHIVE=$SOFTWARE.$VERSION.tar.gz
  URL=http://weblogo.berkeley.edu/release/$ARCHIVE
elif [[ $VERSION == "3.3" ]]
then
  ARCHIVE=$SOFTWARE-$VERSION.tar.gz
  URL=http://weblogo.googlecode.com/files/$ARCHIVE
fi

# If archive was previously downloaded, use the local one, otherwise get it from remote site
if [[ -f ${!INSTALL_HOME}/archive/$ARCHIVE ]]
then
  echo "Archive $ARCHIVE already in ${!INSTALL_HOME}/archive/: using it..."
  cp -a ${!INSTALL_HOME}/archive/$ARCHIVE .
else
  echo "Archive $ARCHIVE not in ${!INSTALL_HOME}/archive/: downloading it..."
  wget $URL
fi
tar zxvf $ARCHIVE

SOFTWARE_DIR=$SOFTWARE-$VERSION

if [[ $VERSION == "2.8.2" ]]
then
  mv $SOFTWARE $SOFTWARE_DIR
  cd $SOFTWARE_DIR
  # Update Perl script shebangs
  sed -i s,"#\!/usr/bin/perl -w,#\!/usr/bin/env perl\\nuse warnings;,g" seqlogo
elif [[ $VERSION == "3.3" ]]
then
  cd $SOFTWARE_DIR
  module load mugqic/python/$PYTHON_VERSION
  python setup.py install
  ln -s weblogo seqlogo
fi

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
prepend-path    PATH                \$root
setenv          WEBLOGO_HOME    \$root
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
