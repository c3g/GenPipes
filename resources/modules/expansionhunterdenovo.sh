#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=ExpansionHunterDenovo
VERSION=0.8.7
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/Illumina/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {

  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR

  mkdir -p build && cd build
  cmake -DCMAKE_BUILD_TYPE=Release ../source
  make

  cd $INSTALL_DOWNLOAD
  mv $SOFTWARE_DIR $INSTALL_DIR/
}


#Module definition to use
module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"HiCPlotter for interaction matrix visualization for Hi-C\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/build
"
}


# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
