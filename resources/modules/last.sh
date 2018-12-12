#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=LAST
VERSION=959
ARCHIVE=${SOFTWARE,,}-$VERSION.zip
ARCHIVE_URL=http://last.cbrc.jp/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE,,}-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE

  cd $SOFTWARE_DIR
  make -j12
  make install prefix=$INSTALL_DIR/$SOFTWARE_DIR
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
prepend-path    PATH                \$root/bin
setenv          LAST_HOME           \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
