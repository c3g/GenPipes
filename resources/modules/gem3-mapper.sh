#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=gem3-mapper
VERSION=3.6.1
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/smarco/${SOFTWARE}/archive/v$VERSION.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  git clone --recursive https://github.com/smarco/gem3-mapper.git -b v${VERSION} $SOFTWARE_DIR
  cd $SOFTWARE_DIR
  ./configure
  make

  # Install software
  cd $INSTALL_DOWNLOAD
  mv $SOFTWARE_DIR $INSTALL_DIR/
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