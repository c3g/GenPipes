#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bio-playground
VERSION=master
ARCHIVE=$SOFTWARE-${VERSION}.zip
ARCHIVE_URL=https://github.com/brentp/$SOFTWARE/archive/$VERSION.zip
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE

  cd $SOFTWARE_DIR/reads-utils
  g++ -O2 -o fastq fastq.cpp

  cd $INSTALL_DOWNLOAD
  mv $SOFTWARE_DIR $INSTALL_DIR/

}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root        $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH        \$root/reads-utils
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

