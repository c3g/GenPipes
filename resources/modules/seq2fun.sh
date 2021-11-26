#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=seq2fun
VERSION=1.2.5
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://www.seq2fun.ca/resources/seq2fun/${SOFTWARE}_v${VERSION}.tar.gz
SOFTWARE_DIR=Seq2Fun-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar -xvzf $ARCHIVE

  cd Seq2Fun/src/
  make clean
  make

  cd $INSTALL_DOWNLOAD
  mv Seq2Fun $INSTALL_DIR/$SOFTWARE_DIR
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
