#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=arriba
VERSION=2.1.0
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/suhrig/${SOFTWARE}/releases/download/v${VERSION}/${SOFTWARE}_v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  mv ${SOFTWARE}_v${VERSION} $SOFTWARE_DIR

  cd $SOFTWARE_DIR
  make -j12

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
prepend-path    PATH                \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
