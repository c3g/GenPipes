#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=cmake
VERSION=3.21.0
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/Kitware/CMake/releases/download/v${VERSION}-rc1/${SOFTWARE}-${VERSION}-rc1.tar.gz
SOFTWARE_DIR=CMake-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd ${SOFTWARE}-${VERSION}-rc1
  ./bootstrap --prefix=$INSTALL_DIR/$SOFTWARE_DIR
  make
  make install
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

