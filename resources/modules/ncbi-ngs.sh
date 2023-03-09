#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=ngs
VERSION=1.3.0
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/ncbi/$SOFTWARE/archive/${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  tar -zxvf $ARCHIVE

  cd $SOFTWARE_DIR 
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR
  make -j12 -C ngs-sdk
  make -j12 -C ngs-java
  make -j12 -C ngs-sdk install
  make -j12 -C ngs-java install
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root ; 
prepend-path    LD_LIBRARY_PATH     \$root/lib64
prepend-path    CLASSPATH           \$root/jar
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
