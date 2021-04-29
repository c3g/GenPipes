#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=sratoolkit
VERSION=2.10.5
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/$VERSION/${SOFTWARE}.${VERSION}-centos_linux64.tar.gz
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  tar -zxvf $ARCHIVE

  mv ${SOFTWARE}.${VERSION}-centos_linux64 $INSTALL_DIR/$SOFTWARE_DIR
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
