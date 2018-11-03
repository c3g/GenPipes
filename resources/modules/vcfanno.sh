#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=vcfanno
VERSION=0.2.9
ARCHIVE=${SOFTWARE}-${VERSION}_linux64	# the archive is actually the binary file
ARCHIVE_URL=https://github.com/brentp/vcfanno/releases/download/v0.2.9/vcfanno_linux64
SOFTWARE_DIR=${SOFTWARE}-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  mkdir -p $SOFTWARE_DIR
  cp $ARCHIVE $SOFTWARE_DIR/${SOFTWARE}-${VERSION}

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
}

#Module definition
module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/
"
}


# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
