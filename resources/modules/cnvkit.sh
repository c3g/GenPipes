#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=CNVkit
VERSION=0.9.5
ARCHIVE=${SOFTWARE,,}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/etal/${SOFTWARE,,}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE,,}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  git clone --recursive https://github.com/etal/${SOFTWARE,,} -b v${VERSION}
  module load mugqic/python/2.7.14
  cd ${SOFTWARE,,}
  pip install -e .

  cd $INSTALL_DOWNLOAD
  mv -i ${SOFTWARE,,} $INSTALL_DIR/$SOFTWARE_DIR
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set              root               $INSTALL_DIR/$SOFTWARE_DIR 
prepend-path     PATH               \$root
prepend-path     PYTHONPATH         \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
