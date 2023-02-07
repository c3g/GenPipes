#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=Quast
VERSION=5.0.2
ARCHIVE=${SOFTWARE,}-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/ablab/${SOFTWARE,}/releases/download/${SOFTWARE,}_${VERSION}/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE,}-$VERSION
PYTHON_VERSION=3.6.5

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  cd $SOFTWARE_DIR
  module load mugqic/python/$PYTHON_VERSION
#  python setup.py install_full --prefix $INSTALL_DIR/$SOFTWARE_DIR
  pip install --upgrade --target=$INSTALL_DIR/$SOFTWARE_DIR $SOFTWARE==$VERSION
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
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

