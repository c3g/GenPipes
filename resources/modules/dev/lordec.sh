#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

SOFTWARE="LoRDEC" 
VERSION="0.5" 
ARCHIVE="LoRDEC-$VERSION.tar.gz"  
ARCHIVE_URL="http://www.atgc-montpellier.fr/download/sources/lordec/LoRDEC-$VERSION.tar.gz" 
SOFTWARE_DIR="$SOFTWARE-$VERSION"  


build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE 
	cd $SOFTWARE_DIR
	module load gcc/4.8.2 BOOST/1.57 
	make install_dep
	make

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
module load gcc/4.8.2 BOOST/1.57 
set             root                $INSTALL_DIR/$SOFTWARE_DIR ;
prepend-path    PATH                \$root ; 

"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
