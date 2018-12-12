#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=VarScan
VERSION=2.4.3
# Replace "." in official version number by "_" in archive version number
ARCHIVE=${SOFTWARE}.v${VERSION}.jar
#ARCHIVE_URL=https://github.com/dkoboldt/${SOFTWARE,,}/releases/download/${VERSION}/$ARCHIVE     # for version < 2.4.3
ARCHIVE_URL=https://github.com/dkoboldt/${SOFTWARE,,}/blob/master/$ARCHIVE                      # for version = 2.4.3
SOFTWARE_DIR=${SOFTWARE}.v${VERSION}

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
setenv          VARSCAN2_HOME        \$root
setenv          VARSCAN2_JAR         \$root/${SOFTWARE}.v${VERSION}.jar
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
