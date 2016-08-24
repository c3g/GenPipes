#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

# BOOST is necessary
module load boost64/1.61.0 
SOFTWARE="NxTrim" 
VERSION="0.4.1"  
ARCHIVE="v$VERSION.tar.gz" 
ARCHIVE_URL="https://github.com/sequencing/NxTrim/archive/$ARCHIVE" 
SOFTWARE_DIR="$SOFTWARE-$VERSION"  

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE 

  cd $SOFTWARE_DIR
	make  

  # Install software
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
module load boost64/1.61.0  ;
set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root ;
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
