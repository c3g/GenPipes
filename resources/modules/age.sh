#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=AGE
VERSION=master
ARCHIVE=${SOFTWARE,,}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/abyzovlab/${SOFTWARE}/archive/${VERSION}.zip
SOFTWARE_DIR=${SOFTWARE}-${VERSION}-20181210

build() {
  cd $INSTALL_DOWNLOAD
  git clone --recursive https://github.com/abyzovlab/AGE.git

  # Install software
  mv ${SOFTWARE} $SOFTWARE_DIR
  cd $SOFTWARE_DIR
  make -j12

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

set              root               $INSTALL_DIR/$SOFTWARE_DIR 
prepend-path     PATH               \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
