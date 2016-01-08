#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=picard
#version 2 or later require JDK1.8
VERSION=2.0.1
ARCHIVE=$SOFTWARE-tools-$VERSION.zip
ARCHIVE_URL=https://github.com/broadinstitute/picard/releases/download/$VERSION/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-tools-$VERSION

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
module-whatis \"$SOFTWARE - require JDK1.8\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          PICARD_HOME         \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
