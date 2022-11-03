#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=gemBS-rs 
VERSION=4.0
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/heathsc/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  git clone --recursive https://github.com/heathsc/gemBS-rs.git -b v${VERSION}

  mv $SOFTWARE $SOFTWARE_DIR
  cd $SOFTWARE_DIR
  make -j12 gemBS_config.mk
  sed -i -e "s|^INSTALL_PREFIX=.*|INSTALL_PREFIX=${INSTALL_DIR}/${SOFTWARE_DIR}|g" gemBS_config.mk
  sed -i -e "s|^GEMBS_INSTALL_ROOT=\$(INSTALL_PREFIX)/.*|GEMBS_INSTALL_ROOT=\$(INSTALL_PREFIX)/|g" gemBS_config.mk
  make -j12 install

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
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
