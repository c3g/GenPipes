#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=HiCPlotter
VERSION=v0.7.3
ARCHIVE=${SOFTWARE}-${VERSION}.zip
ARCHIVE_URL=https://github.com/kcakdemir//$SOFTWARE/archive/master.zip
SOFTWARE_DIR=${SOFTWARE}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {

  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE
  mv ${SOFTWARE}-master $SOFTWARE_DIR

  mv -i $SOFTWARE_DIR $INSTALL_DIR/

  ## add Shebang to make -j12 it easier to call:
  cd $INSTALL_DIR/$SOFTWARE_DIR
  sed -i '1i#!/usr/bin/env python\n' ${SOFTWARE}.py
  chmod 775 ${SOFTWARE}.py

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
prepend-path    PATH                \$root
"
}


# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
