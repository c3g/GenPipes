#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=hmmcopy_utils
COMMIT=29a8d1d
VERSION=master-${COMMIT}
ARCHIVE=${SOFTWARE}-${VERSION}.zip
ARCHIVE_URL=https://github.com/shahcompbio/${SOFTWARE}/archive/refs/heads/master.zip
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd ${INSTALL_DOWNLOAD}

  git clone https://github.com/shahcompbio/${SOFTWARE}.git
  cd ${SOFTWARE}
  git checkout ${COMMIT}
  cmake .
  make

  cd ${INSTALL_DOWNLOAD}
  mv -i ${SOFTWARE} ${INSTALL_DIR}/${SOFTWARE_DIR}
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
