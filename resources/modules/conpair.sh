#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=Conpair
VERSION=0.2.1
ARCHIVE=${SOFTWARE}-${VERSION}.zip
ARCHIVE_URL=https://github.com/mareikejaniak/$SOFTWARE/archive/refs/tags/${VERSION}.zip
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

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
module-whatis \"$SOFTWARE\"

set              root               $INSTALL_DIR/$SOFTWARE_DIR 
prepend-path	 PATH               \$root/scripts
setenv		 CONPAIR_DIR        \$root 
setenv		 CONPAIR_SCRIPTS    \$root/scripts 
setenv		 CONPAIR_DATA       \$root/data 
prepend-path	 PYTHONPATH         \$root/modules
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
