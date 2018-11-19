#!/bin/bash
# Exit immediately on error
set -eu -o pipefail


SOFTWARE=Stacks
VERSION=1.46
ARCHIVE=${SOFTWARE,}-${VERSION}.tar.gz
ARCHIVE_URL=http://catchenlab.life.illinois.edu/${SOFTWARE,}/source/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE,}-${VERSION}

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD

build() {
  cd $INSTALL_DOWNLOAD
  tar xfvz $ARCHIVE

  cd $SOFTWARE_DIR
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR --enable-sparsehash --with-sparsehash-include-path=/cvmfs/soft.mugqic/CentOS6/software/sparsehash_libs/sparsehash-2.0.2/include
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
prepend-path    LD_LIBRARY_PATH     /cvmfs/soft.mugqic/CentOS6/software/sparsehash_libs/sparsehash-2.0.2/lib
prepend-path    CPLUS_INCLUDE_PATH  /cvmfs/soft.mugqic/CentOS6/software/sparsehash_libs/sparsehash-2.0.2/include
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
