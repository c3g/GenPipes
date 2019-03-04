#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=SvABA
VERSION=0.2.2-beta
#VERSION=0.2.1
ARCHIVE=${SOFTWARE,,}-${VERSION}.zip
#ARCHIVE=${SOFTWARE,,}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/walaj/svaba/archive/master.zip
#ARCHIVE_URL=https://github.com/walaj/${SOFTWARE,,}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE,,}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  git clone --recursive https://github.com/walaj/${SOFTWARE,,}
  #git clone --recursive https://github.com/walaj/${SOFTWARE,,} -b ${VERSION}

  # Install
  cd ${SOFTWARE,,}
  ./configure
  make -j12
  make install

  cd $INSTALL_DOWNLOAD
  mv -i ${SOFTWARE,,} $INSTALL_DIR/$SOFTWARE_DIR
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set              root               $INSTALL_DIR/$SOFTWARE_DIR 
prepend-path     PATH               \$root/bin
setenv           SVABA_HOME         \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
