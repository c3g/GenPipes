#!/bin/bash
# Exit immediately on error
set -eu -o pipefail



SOFTWARE="SSPACE-LongRead" 
VERSION="1-1"
ARCHIVE="40$SOFTWARE""_""v$VERSION.tar.gz" 
ARCHIVE_URL="http://www.baseclear.com/base/download/$ARCHIVE"
SOFTWARE_DIR=$SOFTWARE"_"v$VERSION 
# http://www.baseclear.com/base/download/40SSPACE-LongRead_v1-1.tar.gz


build(){ 
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE 

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
setenv          SSPACE_LONGREAD_HOME \$root ; 
prepend-path    PATH                \$root ; 
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
