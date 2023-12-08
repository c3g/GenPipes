#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bs_call
VERSION=2.0.3
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/heathsc/${SOFTWARE}/archive/v$VERSION.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  git clone --recursive https://github.com/heathsc/bs_call.git
  mv $SOFTWARE $SOFTWARE_DIR
  cd $SOFTWARE_DIR
  module load gsl
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR
  make all
  cd ..
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
prepend-path    PATH                \$root/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@