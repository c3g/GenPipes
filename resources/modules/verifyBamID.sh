#!/bin/bash
# Exit immediately on error
set -eu -o pipefail


SOFTWARE=verifyBamID
VERSION=1.1.3
ARCHIVE=${SOFTWARE}_${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/statgen/${SOFTWARE}/releases/download/v${VERSION}/${SOFTWARE}LibStatGen.${VERSION}.tgz
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD

build() {
  cd $INSTALL_DOWNLOAD
  tar -zxvf $ARCHIVE

  cd ${SOFTWARE}_${VERSION}
  make -j12
  make install INSTALLDIR=$INSTALL_DIR/$SOFTWARE_DIR

  cd $INSTALL_DIR/$SOFTWARE_DIR
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
