#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=STAR-Fusion
VERSION=1.10.0
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/${SOFTWARE}/${SOFTWARE}/releases/download/${SOFTWARE}-v${VERSION}/${SOFTWARE}.v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION
MODULE_PERL=mugqic/perl/5.34.0
MODULE_R=mugqic/R_Bioconductor/4.2.1_3.15

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  mv ${SOFTWARE}-v${VERSION} $SOFTWARE_DIR

  module load $MODULE_PERL $MODULE_R
  cd $SOFTWARE_DIR
  make -j12

  # Install software
  cd $INSTALL_DOWNLOAD
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
prepend-path    PATH                \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
