#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

SOFTWARE=qualimap 

# VERSION="v2.0.1"
# ARCHIVE=$SOFTWARE"_"$VERSION.zip
# ARCHIVE_URL=http://qualimap.bioinfo.cipf.es/release/$ARCHIVE
# SOFTWARE_DIR=$SOFTWARE"_"$VERSION


VERSION="2.1.1" 
ARCHIVE=$SOFTWARE"_v"$VERSION.zip  
ARCHIVE_URL="https://bitbucket.org/kokonech/qualimap/downloads/$ARCHIVE"
SOFTWARE_DIR=$SOFTWARE"_v"$VERSION  


# https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.1.1.zip

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE  ## TO BE MODIFIED WITH SPECIFIC COMMAND


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
prepend-path    PATH                \$root ;
setenv          QUALIMAP_JAR        \$root/qualimap.jar ;  
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@


