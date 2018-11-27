#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=freebayes
VERSION=1.2.0
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/ekg/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  git clone --recursive git://github.com/ekg/freebayes.git -b v1.2.0

  cd $SOFTWARE
  make

  # Install software
  cd $INSTALL_DOWNLOAD
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
prepend-path    PATH                \$root/bin/
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
