#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=gridss
VERSION=2.13.2
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/PapenfussLab/${SOFTWARE}/releases/download/v${VERSION}/gridss-${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD

  # Install software
  mkdir -p $SOFTWARE_DIR
  tar zxvf $ARCHIVE
  mv $SOFTWARE_DIR $INSTALL_DIR/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE - require JDK1.8\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          GRIDSS_HOME         \$root
setenv          GRIDSS_JAR          \$root/gridss-2.13.2-gridss-jar-with-dependencies.jar
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

