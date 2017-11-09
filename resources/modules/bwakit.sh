#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bwakit
VERSION=0.7.12
ARCHIVE=${SOFTWARE}-${VERSION}_x64-linux.tar.bz2
ARCHIVE_URL=https://downloads.sourceforge.net/project/bio-bwa/${SOFTWARE}/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar jxvf $ARCHIVE

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i bwa.kit $INSTALL_DIR/$SOFTWARE_DIR
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
