#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=gsort
VERSION=0.1.4
ARCHIVE=${SOFTWARE}-${VERSION}_linux_amd64
ARCHIVE_URL=https://github.com/brentp/${SOFTWARE}/releases/download/v${VERSION}/${SOFTWARE}_linux_amd64
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD

  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR  
  cp $ARCHIVE $INSTALL_DIR/$SOFTWARE_DIR/$SOFTWARE

  chmod 775 $INSTALL_DIR/$SOFTWARE_DIR/$SOFTWARE
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
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
