#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=StringTie
VERSION=1.3.5
ARCHIVE=${SOFTWARE,,}-$VERSION.tar.gz
ARCHIVE_URL=http://ccb.jhu.edu/software/${SOFTWARE,,}/dl/${ARCHIVE}
SOFTWARE_DIR=${SOFTWARE,,}-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  make -j12 release

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
