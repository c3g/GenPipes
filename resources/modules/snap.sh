#!/bin/bash

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# You should probably also delete this commented-out header and the ## comments
################################################################################


#
# Software_name snap.
#

SOFTWARE=snap
VERSION=2013-11-29
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=http://korflab.ucdavis.edu/Software/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # the extraction creates a folder named $SOFTWARE
  # so let's rename it to $SOFTWARE_DIR
  mv $SOFTWARE $SOFTWARE_DIR
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
                      
set             root            $INSTALL_DIR/$SOFTWARE_DIR
repend-path     PATH            \$root/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
