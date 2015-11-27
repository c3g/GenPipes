#!/bin/bash
# Exit immediately on error
set -eu -o pipefail


SOFTWARE=verifyBamID
VERSION=1.1.2  
ARCHIVE=verifyBamIDLibStatGen.1.1.2.tgz
ARCHIVE_URL=https://github.com/statgen/verifyBamID/releases/download/v1.1.2/$ARCHIVE
SOFTWARE_DIR=verifyBamID-1.1.2 

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD

build() {
  cd $INSTALL_DOWNLOAD
  tar -xvzf $ARCHIVE  
  # Transfer immediately from tmp to INSTALL_DIR
  mv -i $INSTALL_DOWNLOAD/$SOFTWARE_DIR $INSTALL_DIR/
  cd  $INSTALL_DIR/$SOFTWARE_DIR
  make install 
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root ;  
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
