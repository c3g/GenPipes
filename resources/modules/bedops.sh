#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bedops
VERSION=v2.4.35
ARCHIVE=${SOFTWARE}_linux_x86_64-${VERSION}.tar.bz2
ARCHIVE_URL=https://github.com/${SOFTWARE}/${SOFTWARE}/releases/download/${VERSION}/${ARCHIVE}
SOFTWARE_DIR=${SOFTWARE}_${VERSION}

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar jxvf $ARCHIVE

  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
  cd $INSTALL_DOWNLOAD
  mv -i bin $INSTALL_DIR/$SOFTWARE_DIR/
}


#Module definition to use
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

