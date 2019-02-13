#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=SVTyper
VERSION=0.7.0
ARCHIVE=${SOFTWARE,,}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/hall-lab/${SOFTWARE,,}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE,,}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  module load mugqic/python/2.7.14
  mkdir -p $SOFTWARE_DIR
  pip install git+https://github.com/hall-lab/svtyper@v${VERSION} --prefix $INSTALL_DIR/$SOFTWARE_DIR
#  git clone --recursive https://github.com/hall-lab/${SOFTWARE,,}.git -b v${VERSION}

  # Install
#  module load mugqic/python/2.7.14
#  cd ${SOFTWARE,,}
#  pip install -e .

  cd $INSTALL_DOWNLOAD
#  mv $SOFTWARE_DIR $INSTALL_DIR/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set              root               $INSTALL_DIR/$SOFTWARE_DIR 
prepend-path     PATH               \$root/bin
prepend-pat      PYTHONPATH         \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
