#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=fastp
VERSION=0.19.7
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/OpenGene/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD

  # What follows is the installation from source (often problematic...)
  git clone --recursive git://github.com/OpenGene/${SOFTWARE}.git -b v$VERSION
  cd $SOFTWARE
  make -j12 
  make install

  # Install software
  cd $INSTALL_DOWNLOAD
  mv $SOFTWARE $INSTALL_DIR/$SOFTWARE_DIR
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          FASTP_HOME          \$root
prepend-path    PATH                \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
