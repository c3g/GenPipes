#!/bin/bash
# Exit immediately on error
set -eu -o pipefail


SOFTWARE=verifyBamID
VERSION=devMaster_20151216
ARCHIVE=${SOFTWARE}_${VERSION}.zip
ARCHIVE_URL=https://github.com/statgen/verifyBamID/archive/master.zip
SOFTWARE_DIR=${SOFTWARE}_${VERSION}

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD

build() {
  echo -e "INSTALL_DOWNLOAD =$INSTALL_DOWNLOAD"
  cd $INSTALL_DOWNLOAD
  unzip ${SOFTWARE}_${VERSION}.zip
  # Transfer immediately from tmp to INSTALL_DIR
  mv -i $INSTALL_DOWNLOAD/verifyBamID-master/ $INSTALL_DIR/$SOFTWARE_DIR/
  cd $INSTALL_DIR/$SOFTWARE_DIR
  mkdir -p bin/
  make cloneLib
  make install INSTALLDIR=$INSTALL_DIR/$SOFTWARE_DIR/bin/ || echo "complete make didn t pass, create module file anyways" 
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
