#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=python
VERSION=3.11.1
# Remove the version last number
LIBVERSION=${VERSION%.[0-9]*}
# Uppercase first P in python
ARCHIVE=${SOFTWARE^}-$VERSION.tgz
ARCHIVE_URL=http://www.python.org/ftp/$SOFTWARE/$VERSION/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE^}-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  # Compile with --enable-unicode=ucs4 to fix error "ImportError: numpy-1.8.1-py2.7-linux-x86_64.egg/numpy/core/multiarray.so: undefined symbol: PyUnicodeUCS2_AsASCIIString"
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR --enable-unicode=ucs4 --enable-loadable-sqlite-extensions --enable-optimizations
  make -j12
  make install

  cd $INSTALL_DIR/$SOFTWARE_DIR/bin
  ln -s python3 python
  ln -s pip3 pip
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          PYTHONHOME          \$root
prepend-path    PATH                \$root/bin
prepend-path    MANPATH             \$root/share/man
prepend-path    LIBRARY_PATH        \$root/lib/
prepend-path    LD_LIBRARY_PATH     \$root/lib/
prepend-path    CPATH               \$root/include:\$root/include/python$LIBVERSION
prepend-path    PYTHONPATH          \$root/lib/python$LIBVERSION/site-packages:\$root/lib/python$LIBVERSION
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
