#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=Manta
VERSION=1.5.0
ARCHIVE=${SOFTWARE,}-${VERSION}.release_src.tar.bz2
ARCHIVE_URL=https://github.com/Illumina/${SOFTWARE,}/releases/download/v${VERSION}/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE,}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  tar -xjf $ARCHIVE

  # Install software
  mkdir build && cd build
  ../${SOFTWARE,}-${VERSION}.release_src/configure --jobs=4 --prefix=$INSTALL_DIR/$SOFTWARE_DIR
  make -j4 install
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set              root               $INSTALL_DIR/$SOFTWARE_DIR 
prepend-path     PATH               \$root/bin
setenv           MANTA_HOME         \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
