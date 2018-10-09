#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=QT
#VERSION=4.8.7
VERSION=5.8.0
ARCHIVE=${SOFTWARE,,}-$VERSION.tar.gz
if [[ $VERSION > 5 ]];
then
  ARCHIVE_URL=http://download.qt.io/official_releases/${SOFTWARE,,}/${VERSION%.*}/${VERSION}/single/${SOFTWARE,,}-everywhere-opensource-src-${VERSION}.tar.gz
else
  ARCHIVE_URL=http://download.qt.io/official_releases/${SOFTWARE,,}/${VERSION%.*}/${VERSION}/${SOFTWARE,,}-everywhere-opensource-src-${VERSION}.tar.gz
fi
SOFTWARE_DIR=${SOFTWARE,,}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd ${SOFTWARE,,}-everywhere-opensource-src-$VERSION
  if [[ $VERSION > 5 ]];
  then
    ./configure -prefix $INSTALL_DIR/$SOFTWARE_DIR
  else
    ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR
  fi
  make -j12
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
