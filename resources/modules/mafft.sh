#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=mafft
VERSION=7.310
ARCHIVE=${SOFTWARE}-$VERSION.tar.gz
ARCHIVE_URL=http://mafft.cbrc.jp/alignment/software/${SOFTWARE}-${VERSION}-with-extensions-src.tgz
SOFTWARE_DIR=${SOFTWARE}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  mv ${SOFTWARE}-${VERSION}-with-extensions $SOFTWARE_DIR

  cd $SOFTWARE_DIR
  OLD_PREFIX="PREFIX = \/usr\/local"
  NEW_PREFIX="PREFIX = $INSTALL_DIR\/$SOFTWARE_DIR"
  cd core
  sed -i "s|$OLD_PREFIX|$NEW_PREFIX|" Makefile
  make -j12
  make install
  cd ../extensions
  sed -i "s|$OLD_PREFIX|$NEW_PREFIX|" Makefile  
  make -j12
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
prepend-path    PATH                \$root/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

