#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=Delly
VERSION=0.7.7
ARCHIVE=${SOFTWARE,,}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/dellytools/${SOFTWARE,,}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE,,}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  git clone --recursive git://github.com/dellytools/delly.git -b v$VERSION

  cd ${SOFTWARE,,}
  make all

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i ${SOFTWARE,,} $INSTALL_DIR/${SOFTWARE_DIR}
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/src  
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
