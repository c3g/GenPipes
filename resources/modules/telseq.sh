#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=TelSeq
VERSION=0.0.2
ARCHIVE=$SOFTWARE-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/zd1/${SOFTWARE,,}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE,,}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar -xzvf $ARCHIVE

  module load mugqic/bamtools/2.5.1

  cd $SOFTWARE_DIR
  cd src
  ./autogen.sh
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR --with-bamtools=$(dirname $(dirname $(which bamtools)))
  make -j12
  make install

#  cd $INSTALL_DOWNLOAD
#  mv $SOFTWARE_DIR $INSTALL_DIR/

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
prepend-path    LD_LIBRARY_PATH     \$root/lib
prepend-path    LIBRARY_PATH        \$root/lib
prepend-path    CPATH               \$root/include
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
