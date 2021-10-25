#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=Nanopolish
VERSION=0.13.3
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/jts/${SOFTWARE,}/archive/refs/tags/v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION
PYTHON_VERSION=3.9.1
PYTHON_SHORT_VERSION=${PYTHON_VERSION:0:3}

build() {
  cd $INSTALL_DOWNLOAD

  git clone --recursive https://github.com/jts/nanopolish.git -b v${VERSION} $SOFTWARE_DIR
  cd $SOFTWARE_DIR
  
  module load mugqic/python/$PYTHON_VERSION
  make

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR/

  # install python dependencies for the Nanoplish helper scripts
  pip install --prefix=$INSTALL_DIR/$SOFTWARE_DIR --ignore-installed -r $INSTALL_DIR/$SOFTWARE_DIR/scripts/requirements.txt 
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
prepend-path    PYTHONPATH          $PYTHONPATH
prepend-path    PYTHONPATH          \$root/lib/python${PYTHON_SHORT_VERSION}
prepend-path    PYTHONPATH          \$root/lib/python${PYTHON_SHORT_VERSION}/site-packages
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
