#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=pandoc
VERSION=1.15.2
ARCHIVE=$SOFTWARE-$VERSION.zip
ARCHIVE_URL=https://s3.amazonaws.com/rstudio-buildtools/$ARCHIVE
#VERSION=1.13.1
#ARCHIVE=$SOFTWARE-$VERSION.zip
#ARCHIVE_URL=https://s3.amazonaws.com/rstudio-buildtools/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE

  # Install software
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
prepend-path    PATH                \$root/linux/debian/x86_64
setenv          LANG                en_US.utf8
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
