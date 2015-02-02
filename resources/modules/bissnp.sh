#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=BisSNP
VERSION=0.82.2
ARCHIVE=$SOFTWARE-$VERSION.jar
ARCHIVE_URL=http://sourceforge.net/projects/bissnp/files/$SOFTWARE-$VERSION/$ARCHIVE/download
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  mkdir $INSTALL_DIR/$SOFTWARE_DIR
  cp $ARCHIVE $INSTALL_DIR/$SOFTWARE_DIR/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          BISSNP_JAR          \$root/$SOFTWARE-$VERSION.jar
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
