#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bwa-mem2
VERSION=2.2.1
ARCHIVE=$SOFTWARE-$VERSION.tar.bz2
ARCHIVE_URL=https://github.com/${SOFTWARE}/${SOFTWARE}/archive/refs/tags/v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION
NOWRAP=1

build() {
  cd $INSTALL_DOWNLOAD
  git clone --recursive https://github.com/bwa-mem2/bwa-mem2 -b v${VERSION} $SOFTWARE_DIR

  cd $SOFTWARE_DIR
  git submodule init
  git submodule update
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
