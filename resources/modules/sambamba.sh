#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=sambamba
VERSION=0.6.8
ARCHIVE=${SOFTWARE}-${VERSION}.gz
ARCHIVE_URL=https://github.com/biod/${SOFTWARE}/releases/download/v${VERSION}/${SOFTWARE}-${VERSION}-linux-static.gz
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  mkdir $SOFTWARE_DIR

  gunzip -c $ARCHIVE > $SOFTWARE_DIR/$SOFTWARE

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
prepend-path    PATH                \$root ; 
setenv          SAMBAMBA_HOME       \$root ;
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
