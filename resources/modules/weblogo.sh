#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

# WebLogo requires a recent version of ghostscript to create PNG and PDF output, and pdf2svg to generate SVG output.
# WebLogo version 3 is written in python. It is necessary to have python 2.5, 2.6 or 2.7 and the extension package numpy installed before WebLogo will run.

SOFTWARE=weblogo
#VERSION=2.8.2
#VERSION=3.3
#VERSION=3.4.1  # not using this becasue missing a shared object file libptf77blas.so.3
VERSION=3.5.0
# If WebLogo version >= 3, specify in which python version it must be intalled
PYTHON_VERSION=2.7.13
if [[ $VERSION == "2.8.2" ]]
then
  ARCHIVE=$SOFTWARE.$VERSION.tar.gz
  ARCHIVE_URL=http://weblogo.berkeley.edu/release/$ARCHIVE
elif [[ $VERSION == "3.3" ]]
then
  ARCHIVE=$SOFTWARE-$VERSION.tar.gz
  ARCHIVE_URL=http://weblogo.googlecode.com/files/$ARCHIVE
elif [[ $VERSION > "3.4" ]]
then
  ARCHIVE=$VERSION-$SOFTWARE.tar.gz
  ARCHIVE_URL=https://github.com/WebLogo/weblogo/archive/$VERSION.tar.gz
fi
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

if [[ $VERSION == "2.8.2" ]]
then
  # For 2.8.2, rename root dir since it is simply 'weblogo'
  mv $SOFTWARE $INSTALL_DIR/$SOFTWARE_DIR
  cd $INSTALL_DIR/$SOFTWARE_DIR
  # Update Perl script shebangs
  sed -i s,"#\!/usr/bin/perl -w,#\!/usr/bin/env perl\\nuse warnings;,g" seqlogo
elif [[ $VERSION == "3.3" ]]
then
  mv $SOFTWARE_DIR $INSTALL_DIR/
  cd $INSTALL_DIR/$SOFTWARE_DIR
  module load mugqic/python/$PYTHON_VERSION
  python setup.py install --prefix $INSTALL_DIR/$SOFTWARE_DIR
  ln -s weblogo seqlogo
elif [[ $VERSION > 3.4 ]]
then
  mv $SOFTWARE_DIR $INSTALL_DIR/
  cd $INSTALL_DIR/$SOFTWARE_DIR
  module load mugqic/python/$PYTHON_VERSION
  pip install --upgrade $SOFTWARE
  ln -s weblogo seqlogo
fi
}

if [[ $VERSION > 3.4 ]]
then
  module_file() {
  echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/../python/Python-${PYTHON_VERSION}
prepend-path    PATH                \$root
prepend-path    PYTHONPATH          \$root/lib/python2.7/site-packages
setenv          WEBLOGO_HOME        \$root
"
  }
else
  module_file() {
  echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
prepend-path    PYTHONPATH          \$root/lib/python2.7/site-packages
setenv          WEBLOGO_HOME        \$root
"
  }
fi

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
