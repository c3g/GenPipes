#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=freebayes
VERSION=1.3.4
#ARCHIVE=${SOFTWARE}-${VERSION}.gz
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
#ARCHIVE_URL=https://github.com/$SOFTWARE/$SOFTWARE/releases/download/v$VERSION/${SOFTWARE}-${VERSION}-linux-static-AMD64.gz
ARCHIVE_URL=https://github.com/$SOFTWARE/$SOFTWARE/releases/download/v$VERSION/${SOFTWARE}-${VERSION}-src.tar.gz
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
#  gunzip -c $ARCHIVE > $SOFTWARE
  git clone --recursive https://github.com/$SOFTWARE/$SOFTWARE.git -b v$VERSION

  cd $SOFTWARE
#  module load mugqic/python/3.7.3 mugqic/htslib/1.11 mugqic/tabix/0.2.6 mugqic/vcflib/1.0.1
  meson build/
  ninja -C build/

  # Install software
  cd $INSTALL_DOWNLOAD
#  mkdir -p $INSTALL_DIR/${SOFTWARE_DIR}
  mv -i $SOFTWARE $INSTALL_DIR/${SOFTWARE_DIR}
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
prepend-path    PATH                \$root/build
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
