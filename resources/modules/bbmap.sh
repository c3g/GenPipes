#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=BBMap
VERSION=38.90
ARCHIVE=$SOFTWARE-$VERSION
ARCHIVE_URL=http://sourceforge.net/projects/bbmap/files/${SOFTWARE}_${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION
NOPATCH=1
NOWRAP=1

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  mv ${SOFTWARE,,} $INSTALL_DIR/$SOFTWARE_DIR
  chmod -R a+rx $INSTALL_DIR/$SOFTWARE_DIR
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
