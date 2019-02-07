#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=Delly
VERSION=0.8.1
ARCHIVE=${SOFTWARE,,}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/dellytools/${SOFTWARE,,}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE,,}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD

  # What follows are intructions to install the binary version of Delly
#  mkdir -p $SOFTWARE_DIR
#  cd $SOFTWARE_DIR
#  wget --no-check-certificate https://github.com/dellytools/${SOFTWARE,,}/releases/download/v${VERSION}/${SOFTWARE,,}_v${VERSION}_linux_x86_64bit
#  wget --no-check-certificate https://github.com/dellytools/${SOFTWARE,,}/releases/download/v${VERSION}/${SOFTWARE,,}_v${VERSION}_parallel_linux_x86_64bit
#  ln -s ${SOFTWARE,,}_v${VERSION}_parallel_linux_x86_64bit ${SOFTWARE,,}
#  chmod 775 ${SOFTWARE,,}*

  # What follows is the installation from source (often problematic...)
  git clone --recursive git://github.com/dellytools/delly.git -b v$VERSION
  mv ${SOFTWARE,,} $SOFTWARE_DIR
  cd $SOFTWARE_DIR
  make all

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
setenv          DELLY_PATH          \$root
prepend-path    PATH                \$root/src
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
