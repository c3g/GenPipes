#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=hmmer
VERSION=3.1b2
#VERSION=3.1b1
#VERSION=2.3.2
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
# Adjust remote download URL according to version first number
if [[ ${VERSION:0:1} == 3 ]]
then
  SUFFIX=3
else
  SUFFIX=""
fi
ARCHIVE_URL=http://selab.janelia.org/software/$SOFTWARE$SUFFIX/$VERSION/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR
  make
  make check
  make install
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
