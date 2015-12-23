#!/bin/bash
# Exit immediately on error
set -eu -o pipefail


SOFTWARE=verifyBamID
VERSION=1.1.2
ARCHIVE=verifyBamID
ARCHIVE_URL=https://github.com/statgen/verifyBamID/releases/download/v1.1.2/verifyBamIDLibStatGen.1.1.2.tgz
SOFTWARE_DIR=verifyBamID_1.1.2

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD

build() {
  cd $INSTALL_DOWNLOAD
  tar -xvzf $ARCHIVE
  # Transfer immediately from tmp to INSTALL_DIR
  mv -i $INSTALL_DOWNLOAD/$SOFTWARE_DIR $INSTALL_DIR/
  cd $INSTALL_DIR/$SOFTWARE_DIR
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR/bin/
  cd verifyBamID
  make install INSTALLDIR=$INSTALL_DIR/$SOFTWARE_DIR/bin/
  cd $INSTALL_DIR/$SOFTWARE_DIR
  make install	|| echo "complete make didn t pass, create module file anyways"
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR/bin
prepend-path    PATH                \$root ;
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
