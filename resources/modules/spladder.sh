#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=spladder
VERSION=1.0.0
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/ratschlab/$SOFTWARE/archive/v$VERSION.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION
SOFTWARE=${SOFTWARE^^[sa]}

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  mv -i $SOFTWARE_DIR $INSTALL_DIR/${SOFTWARE_DIR^^[sa]}

  SOFTWARE_DIR=${SOFTWARE_DIR^^[sa]}
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/python
setenv		SPLADDER_HOME       \$root/python
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
