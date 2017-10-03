#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=juicer
VERSION=0.7.0
ARCHIVE=${SOFTWARE}_tools_${VERSION}.jar
ARCHIVE_URL=http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/${ARCHIVE}
SOFTWARE_DIR=${SOFTWARE}_$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
  cp -i $ARCHIVE $INSTALL_DIR/$SOFTWARE_DIR/  
}


#Module definition to use
module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"Juicer tools for Hi-C\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
setenv          ${SOFTWARE}_JAR     \$root/${ARCHIVE}

"
}


# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
