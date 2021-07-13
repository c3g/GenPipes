#!/bin/bash

SOFTWARE=jellyfish
VERSION=2.3.0
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/gmarcais/${SOFTWARE^}/releases/download/v${VERSION}/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  autoreconf -i
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR
  make -j12
  make install

  # Add permissions and install software
  cd $INSTALL_DOWNLOAD
  chmod -R ug+rwX,o+rX $INSTALL_DIR/$SOFTWARE_DIR
  mv -i $ARCHIVE ${!INSTALL_HOME}/archive
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

