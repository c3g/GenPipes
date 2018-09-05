#!/bin/sh
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=FLASH
VERSION=1.2.11
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://sourceforge.net/projects/flashpage/files/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar xzvf $ARCHIVE

  # Install software
  cd $SOFTWARE_DIR
  make -j12

  # Add permissions and install software
  chmod -R 775 *
  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR
}

# Module file
module_file() {
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE-$VERSION\" ;
}
module-whatis \"$SOFTWARE\" ; 

set     root            $INSTALL_DIR/$SOFTWARE_DIR
setenv  FLASH_HOME      \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

