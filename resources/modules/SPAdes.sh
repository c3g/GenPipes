#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=SPAdes
VERSION=3.13.0
ARCHIVE=${SOFTWARE}-${VERSION}-Linux.tar.gz
#ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=http://cab.spbu.ru/files/release${VERSION}/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE}-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar xzvf $ARCHIVE

#  cd $SOFTWARE_DIR
#  PREFIX=$INSTALL_DIR/$SOFTWARE_DIR ./spades_compile.sh

  mv $SOFTWARE_DIR-Linux $SOFTWARE_DIR

  mv -i $SOFTWARE_DIR $INSTALL_DIR/
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
