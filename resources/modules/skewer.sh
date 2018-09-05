#!/bin/bash
# Exit immediately on error
set -eu -o pipefail


# Have to work from binary because the source code provided in the tarball does not allow to install software anywhere else that /usr/local/bin/
SOFTWARE=skewer
VERSION=0.2.2
ARCHIVE=${SOFTWARE}-${VERSION}-linux-x86_64
ARCHIVE_URL=https://sourceforge.net/projects/${SOFTWARE}/files/Binaries/${ARCHIVE}
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  mkdir -p ${SOFTWARE_DIR}
  cp -p $ARCHIVE ${SOFTWARE_DIR}/skewer
  chmod a+x ${SOFTWARE_DIR}/skewer

  mv ${SOFTWARE_DIR} ${INSTALL_DIR}/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          SKEWER_HOME         \$root
prepend-path    PATH                \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
