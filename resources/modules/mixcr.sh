#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=MiXCR
VERSION=3.0.5
ARCHIVE=${SOFTWARE,,}-${VERSION}.zip
ARCHIVE_URL=https://github.com/milaboratory/${SOFTWARE,,}/releases/download/v${VERSION}/${ARCHIVE}
SOFTWARE_DIR=${SOFTWARE,,}-$VERSION

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
prepend-path    HOME                \$root
setenv          MIXCR_HOME          \$root
setenv          MIXCR_JAR           \$root/${SOFTWARE,,}.jar
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

