#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=openssl
VERSION=1.1.0e
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://www.openssl.org/source/${SOFTWARE}-${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar -zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  ./config --prefix=${INSTALL_DIR}/$SOFTWARE_DIR --openssldir=${INSTALL_DIR}/$SOFTWARE_DIR
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

set             root                            $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                            \$root/bin
prepend-path --delim \" \"        LDFLAGS         -L\$root/lib
prepend-path --delim \" \"        CPPFLAGS        -I\$root/include
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

