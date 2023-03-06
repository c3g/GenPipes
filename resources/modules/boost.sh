#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=boost
VERSION=1.66.0
ARCHIVE=$SOFTWARE-$VERSION.tar.bz2
ARCHIVE_URL=https://sourceforge.net/projects/${SOFTWARE}/files/${SOFTWARE}/${VERSION}/${SOFTWARE}_${VERSION//./_}.tar.bz2
SOFTWARE_DIR=$SOFTWARE-$VERSION

#https://sourceforge.net/projects/boost/files/boost-jam/3.1.18/boost-jam-3.1.18.tgz/download

build() {
  cd $INSTALL_DOWNLOAD
  tar jxvf $ARCHIVE

  mv ${SOFTWARE}_${VERSION//./_} $SOFTWARE_DIR

  # Compile software
  cd $SOFTWARE_DIR
  ./bootstrap.sh --prefix=$INSTALL_DIR/$SOFTWARE_DIR --exec-prefix=$INSTALL_DIR/$SOFTWARE_DIR --libdir=$INSTALL_DIR/$SOFTWARE_DIR/lib --includedir=$INSTALL_DIR/$SOFTWARE_DIR/include --with-python-version=2.7
  ./b2

  # Install software
  ./b2 install

}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"Adds ${SOFTWARE^^} library to your environment\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    LD_LIBRARY_PATH     \$root/lib
prepend-path    LIBRARY_PATH        \$root/lib
prepend-path    CPATH               \$root/include
setenv          BOOST_ROOT          \$root/
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
