#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=FastTree
VERSION=2.1.10
ARCHIVE=${SOFTWARE}-${VERSION}.c
ARCHIVE_URL=http://www.microbesonline.org/${SOFTWARE,,}/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  gcc -O3 -finline-functions -funroll-loops -Wall -o ${SOFTWARE} ${ARCHIVE} -lm
  gcc -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o ${SOFTWARE}MP ${ARCHIVE} -lm

  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
  mv ${SOFTWARE} ${SOFTWARE}MP $INSTALL_DIR/$SOFTWARE_DIR/

  chmod -R 775 $INSTALL_DIR/$SOFTWARE_DIR/*
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
