#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=sparsehash
VERSION=2.0.3
ARCHIVE=${SOFTWARE}-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/$SOFTWARE/$SOFTWARE/archive/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd ${SOFTWARE}-$SOFTWARE_DIR
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

set             root                            $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                            \$root/
prepend-path --delim \" \"      LDFLAGS         -L\$root/lib
prepend-path --delim \" \"      CPPFLAGS        -I\$root/include
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

