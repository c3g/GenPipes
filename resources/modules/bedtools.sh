#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bedtools
VERSION=2.27.0
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
#ARCHIVE_URL=https://bedtools.googlecode.com/files/BEDTools.v${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/arq5x/bedtools2/releases/download/v$VERSION/$ARCHIVE	# for v2.18 and above
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # Rename archive root directory since version is missing
  mv bedtools2 $SOFTWARE_DIR			# for v2.18 and above
#  mv BEDTools-Version-2.16.2 $SOFTWARE_DIR	# for v2.16.2
  cd $SOFTWARE_DIR
  make -j12

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
prepend-path    PATH                \$root/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
