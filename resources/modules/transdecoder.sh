#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=TransDecoder
VERSION=5.7.0
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/${SOFTWARE}/${SOFTWARE}/archive/refs/tags/${SOFTWARE}-v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd ${SOFTWARE}-${SOFTWARE}-v${VERSION}
  make -j12 all

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i ${SOFTWARE}-${SOFTWARE}-v${VERSION} $INSTALL_DIR/${SOFTWARE_DIR}
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - ${SOFTWARE^} \"
}
module-whatis \"${SOFTWARE^}\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
