#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=lofreq
VERSION=2.1.5
ARCHIVE=${SOFTWARE}_star-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/CSB5/${SOFTWARE}/raw/master/dist/${ARCHIVE}
SOFTWARE_DIR=LoFreq-$VERSION
MODULE_HTSLIB=mugqic/htslib/1.14

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  module load $MODULE_HTSLIB

  cd ${SOFTWARE}_star-${VERSION}
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR --with-htslib=$HTSLIB_LIBRARY_DIR
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
prepend-path    PYTHONPATH          \$root/lib/python2.7/site-packages/

"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

