#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=vcflib
VERSION=1.0.0
ARCHIVE=$SOFTWARE-${VERSION%%-*}.tar.gz
ARCHIVE_URL=https://github.com/$SOFTWARE/$SOFTWARE/archive/v$VERSION.tar.gz
SOFTWARE_DIR=$SOFTWARE-${VERSION%%-*}

# Here we chop the VERSION to only keep relevant digits
# before it is passed to install_module.sh
VERSION=${VERSION%%-*}

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  git clone --recursive git://github.com/ekg/vcflib.git $SOFTWARE_DIR

  # Install software
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

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
