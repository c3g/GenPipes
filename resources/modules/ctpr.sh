#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=CTPR
VERSION=master_20190115
ARCHIVE=$SOFTWARE-${VERSION}.zip
ARCHIVE_URL=https://github.com/wonilchung/$SOFTWARE/archive/master.zip
SOFTWARE_DIR=$SOFTWARE-${VERSION}

build() {
  cd $INSTALL_DIR
  git clone https://github.com/wonilchung/$SOFTWARE
  mv $SOFTWARE $SOFTWARE_DIR
  chmod a+x $SOFTWARE_DIR/*
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
