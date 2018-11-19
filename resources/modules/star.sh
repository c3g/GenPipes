#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=star
VERSION=2.5.4b
ARCHIVE=$VERSION.tar.gz # for 2.5.0b and newer
#ARCHIVE=${SOFTWARE^^}_$VERSION.tar.gz # for 2.5.0a and older 
ARCHIVE_URL=https://github.com/alexdobin/STAR/archive/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE^^}_$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  
  # Remove "STAR-" prefix from top directory name
#  mv ${SOFTWARE^^}-$SOFTWARE_DIR $SOFTWARE_DIR # for 2.5.0a and older
  cd ${SOFTWARE^^}-$VERSION/source
  make -j12 STAR STARlong

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i ${SOFTWARE^^}-$VERSION $INSTALL_DIR/$SOFTWARE_DIR
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/source
module load muqgic/gcc/4.9.3
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
