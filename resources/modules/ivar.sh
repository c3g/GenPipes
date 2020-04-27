#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=ivar
VERSION=1.2.1
ARCHIVE=$SOFTWARE-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/andersen-lab/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar -xzvf $ARCHIVE

  cd $SOFTWARE_DIR
  module load mugqic/htslib/1.9
  ./autogen.sh
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR
  make
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
module load mugqic/samtools/1.9 mugqic/htslib/1.9
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
