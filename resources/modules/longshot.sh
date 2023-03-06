#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=longshot
VERSION=0.4.1
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/pjedge/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar -xzvf $ARCHIVE

  cd $SOFTWARE_DIR
  C_INCLUDE_PATH=$INSTALL_DIR/$SOFTWARE_DIR/include LIBRARY_PATH=$INSTALL_DIR/$SOFTWARE_DIR/lib cargo install --verbose --path . --root $INSTALL_DIR/$SOFTWARE_DIR
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
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

