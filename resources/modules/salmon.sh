#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=salmon
VERSION=1.3.0
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
# Binary
#ARCHIVE_URL=https://github.com/COMBINE-lab/${SOFTWARE}/releases/download/v${VERSION}/${SOFTWARE}-${VERSION}_linux_x86_64.tar.gz
# Source code
ARCHIVE_URL=https://github.com/COMBINE-lab/${SOFTWARE}/archive/refs/tags/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE^}-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
#  # If we get the binary
#  mv ${SOFTWARE}-${VERSION}_linux_x86_64 $INSTALL_DIR/$SOFTWARE_DIR

  # If we compile from source
  cd ${SOFTWARE_DIR,}
  # For salmon 1.3.0 add the following sed command
  sed -i -e 's/namespace salmon/#include <string>\n\nnamespace salmon/' include/BAMUtils.hpp
  mkdir -p build
  cd build
  rm -rf *
  cmake \
   -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/$SOFTWARE_DIR \
   -DNO_IPO=TRUE \
   ..
  make -j12
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
prepend-path    LD_LIBRARY_PATH     \$root/lib
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
