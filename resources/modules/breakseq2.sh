#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=breakseq2
VERSION=2.2
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/bioinform/$SOFTWARE/archive/$VERSION.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION
PYTHON_MODULE=mugqic/python/2.7.14

build() {
  cd $INSTALL_DOWNLOAD

  module load $PYTHON_MODULE
  cd $INSTALL_DOWNLOAD
  pip install --upgrade --target=$INSTALL_DIR/$SOFTWARE_DIR $ARCHIVE_URL
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
prepend-path    PYTHONPATH          \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
