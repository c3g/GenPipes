#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=mysql
VERSION=8.0.32
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/${SOFTWARE}/${SOFTWARE}-server/archive/refs/tags/${ARCHIVE}
SOFTWARE_DIR=${SOFTWARE}-${VERSION}
MODULE_BOOST=mugqic/boost/1.77.0
build() {
  cd $INSTALL_DOWNLOAD
  tar -xvf ${ARCHIVE}
  cd ${SOFTWARE}-server-${SOFTWARE_DIR}
  rm -rf bld
  mkdir bld
  cd bld
  module load ${MODULE_BOOST}
  cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}/${SOFTWARE_DIR} ..
  make -j12
  make install
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE - require JDK1.8\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@