#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=kallisto 
VERSION=0.48.0
ARCHIVE=v${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/pachterlab/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE}-${VERSION}  

build() {
  cd ${INSTALL_DOWNLOAD}
  rm -rf ${SOFTWARE}
  git clone https://github.com/pachterlab/${SOFTWARE}.git -b v${VERSION}

  cd ${SOFTWARE}

  mkdir build
  cd build
  cmake -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}/${SOFTWARE_DIR} ..
  make
  make install
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - ${SOFTWARE} \"
}
module-whatis \"${SOFTWARE}\"

set             root                ${INSTALL_DIR}/${SOFTWARE_DIR}
prepend-path    PATH                \$root/bin ; 
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
