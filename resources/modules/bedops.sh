#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bedops
VERSION=v2.4.40
ARCHIVE=${SOFTWARE}_linux_x86_64-${VERSION}.tar.bz2
ARCHIVE_URL=https://github.com/$SOFTWARE/$SOFTWARE/archive/refs/tags/$VERSION.tar.gz
SOFTWARE_DIR=${SOFTWARE}_${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  git clone https://github.com/bedops/bedops.git -b $VERSION $SOFTWARE_DIR

  cd $SOFTWARE_DIR
  make all -j12
  make install_all

  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
}


#Module definition to use
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

