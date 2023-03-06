#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=medaka
VERSION=1.0.3
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/nanoporetech/${SOFTWARE}/archive/v$VERSION.tar.gz
SOFTWARE_DIR=${SOFTWARE}-$VERSION
PYTHON_VERSION=3.8.5
NOPATCH=1

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  cd $SOFTWARE_DIR
  module load mugqic/python/$PYTHON_VERSION
#  python setup.py install --prefix $INSTALL_DIR/$SOFTWARE_DIR
  pip install --no-index --upgrade --target=$INSTALL_DIR/$SOFTWARE_DIR $SOFTWARE==$VERSION
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - ${SOFTWARE^} \"
}
module-whatis \"${SOFTWARE^}\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
prepend-path    PYTHONPATH          $PYTHONPATH
prepend-path    PYTHONPATH          \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

