#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

SOFTWARE="pbsuite"  
VERSION="15.2.20"  
ARCHIVE="PBSuite_$VERSION.p1.tgz"
ARCHIVE_URL="http://downloads.sourceforge.net/project/pb-jelly/$ARCHIVE" 
SOFTWARE_DIR="PBSuite_$VERSION"


# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar -xvf $ARCHIVE  

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

# Dependencies: python, networkx, blasr
prereq mugqic/python/2.7.8
prereq mugqic/smrtanalysis/2.3.0.140936.p2

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          SWEETPATH           \$root;
prepend-path    PATH                \$root/bin; 
prepend-path    PYTHONPATH          \$root; 

"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
