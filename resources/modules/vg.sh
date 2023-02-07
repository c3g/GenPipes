#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=vg
VERSION=1.28.0
ARCHIVE=${SOFTWARE}-v${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/vgteam/${SOFTWARE}/releases/download/v${VERSION}/${SOFTWARE}-v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  wget https://github.com/vgteam/${SOFTWARE}/releases/download/v${VERSION}/${SOFTWARE}
  
  mkdir -p $SOFTWARE_DIR
  mv ${SOFTWARE} ${SOFTWARE_DIR}/
  mv $SOFTWARE_DIR ${INSTALL_DIR}/
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

