#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=MUSCLE
VERSION=3.8.31
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=http://www.drive5.com/${SOFTWARE,,}/downloads${VERSION}/${SOFTWARE,,}${VERSION}_i86linux64.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # Install software
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR/bin
  mv -i muscle${VERSION}_i86linux64 $INSTALL_DIR/$SOFTWARE_DIR/bin/
  cd $INSTALL_DIR/$SOFTWARE_DIR/bin/
  ln -s muscle${VERSION}_i86linux64 muscle
  chmod -R 775 $INSTALL_DIR/$SOFTWARE_DIR
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
