#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=emboss
VERSION=6.4.0
ARCHIVE=${SOFTWARE^^}-$VERSION.tar.gz
ARCHIVE_URL=ftp://emboss.open-bio.org/pub/${SOFTWARE^^}/old/$VERSION/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE^^}-$VERSION

# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR
  make
  make install

#  # Install software
#  cd $INSTALL_DOWNLOAD
#  mv -i $SOFTWARE_DIR $INSTALL_DIR/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
<<<<<<< HEAD
prepend-path    PATH                \$root/bin 
=======
prepend-path    PATH                \$root/bin ;  ## TO BE ADDED IF NECESSARY
prepend-path    PATH                \$root/other_tools/bin ;  ## TO BE ADDED AND MODIFIED IF NECESSARY
setenv          ${SOFTWARE}_JAR     \$root/$SOFTWARE-$VERSION.jar ;  ## TO BE ADDED AND MODIFIED IF NECESSARY
>>>>>>> 8fd6a71dd42f95b41106e0e7ff5f1b6e00066d63
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
