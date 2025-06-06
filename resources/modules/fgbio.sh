#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=fgbio
VERSION=2.1.0
ARCHIVE=${SOFTWARE}-${VERSION}.jar
ARCHIVE_URL=https://github.com/fulcrumgenomics/${SOFTWARE}/releases/download/${VERSION}/${ARCHIVE}
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD

  mkdir -p $SOFTWARE_DIR
  cp $ARCHIVE $SOFTWARE_DIR/

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
setenv          FGBIO_JAR           \$root/$ARCHIVE
prepend-path    PATH                \$root/
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
