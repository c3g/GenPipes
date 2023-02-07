#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=variantbam
SOFTREAL=VariantBam
VERSION=1.4.3
ARCHIVE=${SOFTREAL}-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/walaj/${SOFTREAL}/archive/refs/tags/${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTREAL}-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  rm -rf ${SOFTREAL}

  git clone --recursive https://github.com/jwalabroad/${SOFTREAL}.git -b $VERSION
  cd ${SOFTREAL}
  ./configure --prefix=$INSTALL_DIR/${SOFTWARE_DIR}
  make -j12
  
  ## set hts to build with libcurl links and hfile_libcurl.c
  cd SeqLib/htslib
  ./configure --enable-libcurl --prefix=$INSTALL_DIR/${SOFTWARE_DIR} 
  ## compile seqlib with libcurl support
  cd ../../ # back to VariantBam main directory
  ./configure LDFLAGS="-lcurl -lcrypto" --prefix=$INSTALL_DIR/${SOFTWARE_DIR}
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
prepend-path    PATH                \$root/bin"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
