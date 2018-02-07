#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=InterOp
VERSION=1.1.1
ARCHIVE=${SOFTWARE}-${VERSION}-Linux-GNU.tar.gz
ARCHIVE_URL=https://github.com/Illumina/${SOFTWARE,,}/releases/download/v${VERSION}/${ARCHIVE}

SOFTWARE_DIR=${SOFTWARE}-$VERSION

# Specific commands to extract�| and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # Install software
  mv -i $SOFTWARE_DIR-Linux-GNU $INSTALL_DIR/$SOFTWARE_DIR
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
