#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=dragmap
VERSION=1.2.1
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/Illumina/${SOFTWARE^^}/archive/refs/tags/${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE^^}-$VERSION

build() {
  cd $INSTALL_DOWNLOAD

  git clone https://github.com/Illumina/DRAGMAP.git -b ${VERSION}
  cd ${SOFTWARE^^}
  make -j12

  # Install software
  cd $INSTALL_DOWNLOAD
  mv ${SOFTWARE^^} $INSTALL_DIR/${SOFTWARE^^}-$VERSION
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/build/release
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
