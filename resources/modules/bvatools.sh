#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bvatools
VERSION=1.6
ARCHIVE=$SOFTWARE-$VERSION.zip
ARCHIVE_URL=https://bitbucket.org/mugqic/$SOFTWARE/downloads/$ARCHIVE
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
setenv          BVATOOLS_HOME       \$root
setenv          BVATOOLS_JAR        \$root/$SOFTWARE-$VERSION-full.jar
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
