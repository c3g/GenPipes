#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=anaconda
VERSION=2-4.0.0
ARCHIVE=${SOFTWARE^}${VERSION}-Linux-x86_64.sh
ARCHIVE_URL=https://repo.continuum.io/archive/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE^}${VERSION%-[0-9]*.[0-9]*.[0-9]*}

# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  bash $ARCHIVE -b -f -p $INSTALL_DIR/$SOFTWARE_DIR
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
