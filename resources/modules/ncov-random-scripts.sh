#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=ncov-random-scripts
VERSION=20210415
ARCHIVE=$SOFTWARE-$VERSION.zip
ARCHIVE_URL=https://github.com/jts/${SOFTWARE}/archive/refs/heads/master.zip
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD

  git clone https://github.com/jts/ncov-random-scripts.git $SOFTWARE_DIR

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

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
if { [ module-info mode load ] } {
    puts        stderr              \"Ncov random scripts need python3 with pysam and parasail...\"
    puts        stderr              \"e.g. loadig mugqic/python/3.7.3 would set the proper environment\"
}
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
