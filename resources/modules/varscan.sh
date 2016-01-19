#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=VarScan
VERSION=2.3.9
# Replace "." in official version number by "_" in archive version number
ARCHIVE=${SOFTWARE}.v${VERSION}.jar
ARCHIVE_URL=http://sourceforge.net/projects/varscan/files/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE}.v${VERSION}

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD

  # Install software
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
  cp -i $ARCHIVE $INSTALL_DIR/$SOFTWARE_DIR/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          VARSCAN_HOME        \$root
setenv          VARSCAN_JAR         \$root/${SOFTWARE}.v${VERSION}.jar
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
