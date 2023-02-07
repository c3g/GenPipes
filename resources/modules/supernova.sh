#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=supernova
VERSION=2.1.1
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
# supernova archive has to be manually downloaded from https://support.10xgenomics.com/de-novo-assembly/software/downloads/latest
# and then stored in $MUGQIC_INSTALL_HOME/archive/ or/and $MUGQIC_INSTALL_HOME_DEV/archive/
ARCHIVE_URL=
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # Move software
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
