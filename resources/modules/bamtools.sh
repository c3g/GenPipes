#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bamtools
VERSION=2.5.1
ARCHIVE=$SOFTWARE-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/pezmaster31/$SOFTWARE/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar -xzvf $ARCHIVE

  cd $SOFTWARE_DIR
  mkdir -p build
  cd build
  cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/$SOFTWARE_DIR ..
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
prepend-path    LIBRARY_PATH        \$root/lib
prepend-path    CPATH               \$root/include
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
