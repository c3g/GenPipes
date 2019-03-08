#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=CMake
VERSION=3.14.0
ARCHIVE=${SOFTWARE,,}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/Kitware/${SOFTWARE}/releases/download/v${VERSION}-rc2/${SOFTWARE,,}-${VERSION}-rc2.tar.gz
SOFTWARE_DIR=${SOFTWARE,,}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  mv ${SOFTWARE,,}-${VERSION}-rc2 $SOFTWARE_DIR

  cd $SOFTWARE_DIR
  ./bootstrap --prefix=$INSTALL_DIR/$SOFTWARE_DIR
  make
  make install

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          DELLY_PATH          \$root
prepend-path    PATH                \$root/src
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

