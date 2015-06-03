#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

SOFTWARE=usearch 
VERSION=6.1.544
ARCHIVE=$SOFTWARE-$VERSION.tar.gz  
ARCHIVE_URL=http://www.drive5.com/usearch/download.html
SOFTWARE_DIR=$SOFTWARE-$VERSION  

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {

  # Download usearch v6.1.544 on http://www.drive5.com/usearch/download.html.
  # Rename the file to 'usearch61' and add the right for execution.
  # Move it into the folder software/
  
  cd $INSTALL_DIR
  mkdir $SOFTWARE_DIR
  mv $MUGQIC_INSTALL_HOME/software/usearch61 $INSTALL_DIR/$SOFTWARE_DIR  ## TO BE ADDED AND MODIFIED IF NECESSARY
  export PATH=$INSTALL_DIR/$SOFTWARE_DIR:$PATH
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          USEARCH61_HOME         \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
