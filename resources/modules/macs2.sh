#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=MACS2
VERSION=2.1.1.20160309
#VERSION=2.1.0.20151222
#VERSION=2.1.0.20140616
PYTHON_VERSION=2.7.11
#PYTHON_VERSION=2.7.8

ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://pypi.python.org/packages/source/M/MACS2/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
# recent version of GCC is needed
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  module load mugqic/python/$PYTHON_VERSION
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR/lib/python2.7/site-packages
  export PYTHONPATH=${PYTHONPATH}:$INSTALL_DIR/$SOFTWARE_DIR/lib/python2.7/site-packages
  python setup.py install --prefix $INSTALL_DIR/$SOFTWARE_DIR

  # Update Python script shebang
  sed -i s,"#\!/.*python,#\!/usr/bin/env python,g" $INSTALL_DIR/$SOFTWARE_DIR/bin/macs2
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
prepend-path    PYTHONPATH          \$root/lib/python2.7/site-packages
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
