#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=ViennaRNA
VERSION=2.4.14
VERSION_FAMILY=${VERSION:0:4}x
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://www.tbi.univie.ac.at/RNA/download/sourcecode/${VERSION_FAMILY//./_}/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR
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
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

