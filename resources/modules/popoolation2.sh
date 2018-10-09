#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=popoolation2
VERSION=1201
ARCHIVE=${SOFTWARE}_${VERSION}.zip
ARCHIVE_URL=https://sourceforge.net/projects/popoolation2/files/${ARCHIVE}
SOFTWARE_DIR=${SOFTWARE}_${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE

  mkdir -p ${INSTALL_DIR}/${SOFTWARE_DIR}
  mv ${SOFTWARE_DIR}/* ${INSTALL_DIR}/${SOFTWARE_DIR}/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                    $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                    \$root/
setenv          POPOOLATION2_HOME       \$root/
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
