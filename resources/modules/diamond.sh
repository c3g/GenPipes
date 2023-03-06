#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=diamond
VERSION=2.1.4
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=http://github.com/bbuchfink/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE}-${VERSION}
MODULE_BLAST=mugqic/blast/2.13.0+

build() {
  cd $INSTALL_DOWNLOAD

  tar zxvf $ARCHIVE
  cd ${SOFTWARE_DIR}

  mkdir -p bin
  cd bin

  module load $MODULE_BLAST
  cmake -DWITH_ZSTD=ON \
    -DBLAST_INCLUDE_DIR=$BLAST_HOME/include/ncbi-tools++ \
    -DBLAST_LIBRARY_DIR=$BLAST_HOME/lib \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}/${SOFTWARE_DIR} \
    ..
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
