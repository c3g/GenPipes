#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

VERSION=0.5.2  
ARCHIVE=$SOFTWARE-$VERSION.tar.gz 
ARCHIVE_URL=http://sourceforge.net/projects/${SOFTWARE}/files/${ARCHIVE} 
SOFTWARE_DIR=$SOFTWARE-$VERSION 


build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE  

  cd $SOFTWARE_DIR
  make -j12  
  

  # Install software
  cd $INSTALL_DOWNLOAD  ## TO BE ADDED AND MODIFIED IF NECESSARY
  mv -i $SOFTWARE_DIR $INSTALL_DIR/  ## TO BE ADDED AND MODIFIED IF NECESSARY
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
prepend-path    LD_LIBRARY_PATH     \$root/bamtools-2.3.0/lib/
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
