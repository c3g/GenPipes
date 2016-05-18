#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

SOFTWARE=bcbio.variation
VERSION=0.2.6
ARCHIVE=$SOFTWARE-$VERSION-standalone.jar 
ARCHIVE_URL=https://github.com/chapmanb/bcbio.variation/releases/download/v${VERSION}/$ARCHIVE
SOFTWARE_DIR=$VERSION 

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD

  mkdir -p $INSTALL_DIR/$VERSION

  mv -i $INSTALL_DOWNLOAD/* $INSTALL_DIR/$VERSION
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                    $INSTALL_DIR/$VERSION
prepend-path    PATH                    \$root ;  
setenv          BCBIO_VARIATION_JAR     \$root/${SOFTWARE}-${VERSION}-standalone.jar ; 
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
