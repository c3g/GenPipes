#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bcbio.variation.recall
VERSION=0.2.6
ARCHIVE=${SOFTWARE//./-}-$VERSION
ARCHIVE_URL=https://github.com/chapmanb/$SOFTWARE/releases/download/v$VERSION/${SOFTWARE//./-}
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the softwa
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  mkdir $SOFTWARE_DIR

  # Install software
  cp $ARCHIVE $SOFTWARE_DIR/$SOFTWARE
  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                            $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                            \$root  
setenv          BCBIO_VARIATION_RECALL_HOME     \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

