#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=vawk
VERSION=0.0.2
#ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE=${SOFTWARE}-${VERSION}.tar.bz2
#ARCHIVE_URL=https://github.com/cc2qe/${SOFTWARE}/archive/${VERSION}.tar.gz
ARCHIVE_URL=https://anaconda.org/bioconda/${SOFTWARE}/${VERSION}/download/noarch/${SOFTWARE}-${VERSION}-py_4.tar.bz2
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
#  tar zxvf $ARCHIVE
  tar jxvf $ARCHIVE

#  mv -i $SOFTWARE_DIR $INSTALL_DIR/
   mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
   cp python-scripts/$SOFTWARE $INSTALL_DIR/$SOFTWARE_DIR/
   chmod 775 $INSTALL_DIR/$SOFTWARE_DIR/$SOFTWARE
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
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
