#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=java
VERSION=openjdk-jdk-17
ARCHIVE=jdk-17_linux-x64_bin.tar.gz
ARCHIVE_URL=https://download.oracle.com/java/17/latest/$ARCHIVE
SOFTWARE_DIR=$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # Install software
  mv -i ${SOFTWARE_DIR/openjdk-/}.0.4.1 $INSTALL_DIR/$SOFTWARE_DIR
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
