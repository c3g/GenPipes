#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=picard
#version 2 or later require JDK1.8
VERSION=2.27.4
ARCHIVE=${SOFTWARE}-${VERSION}.jar
ARCHIVE_URL=https://github.com/broadinstitute/picard/releases/download/$VERSION/${SOFTWARE}.jar
SOFTWARE_DIR=$SOFTWARE-$VERSION
MODULE_JAVA=mugqic/java/openjdk-jdk-19.0.1

build() {
  cd $INSTALL_DOWNLOAD

  # Install software
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
  cp $ARCHIVE $INSTALL_DIR/$SOFTWARE_DIR/${SOFTWARE}.jar
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE - require JDK1.8\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          PICARD_HOME         \$root/build/libs
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

