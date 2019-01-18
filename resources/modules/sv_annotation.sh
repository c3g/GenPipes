#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=simple_sv_annotation
VERSION=1.0.0
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/AstraZeneca-NGS/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # Install software
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set              root               $INSTALL_DIR/$SOFTWARE_DIR 
setenv           SVANNOT_PATH       \$root
prepend-path     PATH               \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
