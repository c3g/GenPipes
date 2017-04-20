#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=python
VERSION=2.7.13
SETUPTOOLS_VERSION=20.9.0
# Remove the version last number
LIBVERSION=${VERSION%.[0-9]*}
# Uppercase first P in python
ARCHIVE=${SOFTWARE^}-$VERSION.tgz
ARCHIVE_URL=http://www.python.org/ftp/$SOFTWARE/$VERSION/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE^}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  # Compile with --enable-unicode=ucs4 to fix error "ImportError: numpy-1.8.1-py2.7-linux-x86_64.egg/numpy/core/multiarray.so: undefined symbol: PyUnicodeUCS2_AsASCIIString"
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR --enable-unicode=ucs4 --with-zlib-dir=/usr/lib64 --with-ensurepip=install
  make -j8
  make install

  echo "General Python installation done.... processing setuptools"

  # Install setuptools => easy_install
  cd $INSTALL_DOWNLOAD
  SETUPTOOLS_ARCHIVE=setuptools-${SETUPTOOLS_VERSION}.tar.gz
  download_archive https://pypi.python.org/packages/source/s/setuptools $SETUPTOOLS_ARCHIVE
  tar zxvf $SETUPTOOLS_ARCHIVE
  cd ${SETUPTOOLS_ARCHIVE/.tar.gz/}
  #SETUPTOOLS_ARCHIVE=v${SETUPTOOLS_VERSION}.tar.gz
  #download_archive https://github.com/pypa/setuptools/archive $SETUPTOOLS_ARCHIVE
  #tar zxvf $SETUPTOOLS_ARCHIVE
  #cd setuptools-$SETUPTOOLS_VERSION
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python setup.py build
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python setup.py install

  EASY_INSTALL_PATH=$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install
  
  # pip
  ${EASY_INSTALL_PATH} pip
  PIP_PATH=$INSTALL_DIR/$SOFTWARE_DIR/bin/pip
  
  #Add permissions
  chmod -R ug+rwX,o+rX-w $INSTALL_DIR/$SOFTWARE_DIR

}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          PYTHON_HOME         \$root
prepend-path    PATH                \$root/bin
prepend-path    MANPATH             \$root/share/man
prepend-path    LIBRARY_PATH        \$root/lib/
prepend-path    LD_LIBRARY_PATH     \$root/lib/
prepend-path    CPATH               \$root/include:\$root/include/python$LIBVERSION
prepend-path    PYTHONPATH          \$root/lib/python$LIBVERSION/site-packages:\$root/lib/python$LIBVERSION
setenv          QIIME_HOME          \$root/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

