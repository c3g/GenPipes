#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=GangSTR
VERSION=2.4.4
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/gymreklab/${SOFTWARE}/releases/download/v${VERSION}/${ARCHIVE}
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
  ./install-gangstr.sh $INSTALL_DIR/$SOFTWARE_DIR

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
