#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=LUMPY-SV
VERSION=0.2.13
ARCHIVE=${SOFTWARE,,}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/arq5x/${SOFTWARE,,}/releases/download/${VERSION}/${SOFTWARE,,}-v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE,,}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  mv ${SOFTWARE,,}-v${VERSION} $SOFTWARE_DIR

  # Install software
  cd $SOFTWARE_DIR
  make -j12

  # Correct lumpyexpress.config
  cd scripts
  sed -i 's/LUMPY_HOME=~\/lumpy-sv/#LUMPY_HOME=~\/lumpy-sv/g' lumpyexpress.config

  cd $INSTALL_DOWNLOAD
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
prepend-path     PATH               \$root/bin
setenv           LUMPY_HOME         \$root
setenv           LUMPY_SCRIPTS      \$root/scripts
prereq mugqic/python/2.7.14 mugqic/samblaster/0.1.24 mugqic/sambamba/0.6.6 mugqic/samtools/1.4.1
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
