#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=python
VERSION=2.7.11
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
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR --enable-unicode=ucs4
  make -j8
  make install

  echo "General Python installation done.... processing setuptools"

  # Install setuptools => easy_install
  cd $INSTALL_DOWNLOAD
  SETUPTOOLS_ARCHIVE=setuptools-${SETUPTOOLS_VERSION}.tar.gz
  download_archive https://pypi.python.org/packages/source/s/setuptools $SETUPTOOLS_ARCHIVE 
  tar zxvf $SETUPTOOLS_ARCHIVE
  cd ${SETUPTOOLS_ARCHIVE/.tar.gz/}
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python setup.py build
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python setup.py install

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
"
}

create_dir() {
  DIR=$1

  # Create directory with permissions if necessary
  if [[ ! -d $DIR ]]
  then
    mkdir -p $DIR
    chmod ug+rwX,o+rX-w $DIR
  fi
}

download_archive() {
  INSTALL_DOWNLOAD=$INSTALL_DIR/tmp
  mkdir -p $INSTALL_DOWNLOAD

    if [[ "$#" -eq 2 ]]
    then
        ARCHIVE_TMP=$2
        ARCHIVE_URL_PREFIX_TMP=$1
        ARCHIVE_URL_TMP=${ARCHIVE_URL_PREFIX_TMP}/${ARCHIVE_TMP}
    else
        ARCHIVE_TMP=$ARCHIVE
        ARCHIVE_URL_TMP=$ARCHIVE_URL
    fi

  # If archive was previously downloaded, use the local one, otherwise get it from remote site
  if [[ -f $ARCHIVE_DIR/$ARCHIVE_TMP ]]
  then
    echo "Archive $ARCHIVE_TMP already in $ARCHIVE_DIR/: using it..."
    cp -a $ARCHIVE_DIR/$ARCHIVE_TMP $INSTALL_DOWNLOAD/
  else
    echo "Archive $ARCHIVE_TMP not in $ARCHIVE_DIR/: downloading it..."
    wget --no-check-certificate $ARCHIVE_URL_TMP --output-document=$INSTALL_DOWNLOAD/$ARCHIVE_TMP
  fi
}

store_archive() {
  ARCHIVE_TMP=$1

  # Store archive if not already present
  if [[ ! -f $ARCHIVE_DIR/$ARCHIVE_TMP ]]
  then
    chmod -R ug+rwX,o+rX-w $INSTALL_DOWNLOAD/$ARCHIVE_TMP
    create_dir $ARCHIVE_DIR
    mv $INSTALL_DOWNLOAD/$ARCHIVE_TMP $ARCHIVE_DIR/
  fi
}

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
if [[ ${1:-} == MUGQIC_INSTALL_HOME ]]
then
  INSTALL_HOME=MUGQIC_INSTALL_HOME
  MODULE_HOME=muqgic
elif [[ ${1:-} == MUGQIC_INSTALL_HOME_ED ]]
then
  INSTALL_HOME=MUGQIC_INSTALL_HOME_ED
  MODULE_HOME=muqgic_ed
else
  INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV
  MODULE_HOME=muqgic_dev
fi

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE
ARCHIVE_DIR=${!INSTALL_HOME}/archive

# Set module directory path by lowercasing $INSTALL_HOME and removing '_install_home' in it
MODULE_DIR=${!INSTALL_HOME}/modulefiles/`echo ${INSTALL_HOME,,} | sed 's/_install_home//'`/$SOFTWARE

echo "Installing $SOFTWARE version ${VERSION}_qiime in \$$INSTALL_HOME..."
echo

# Abort if software and/or module are already installed
if [[ -e $INSTALL_DIR/$SOFTWARE_DIR || -e $MODULE_DIR/${VERSION}_qiime ]]
then
  echo "$INSTALL_DIR/$SOFTWARE_DIR and/or $MODULE_DIR/${VERSION}_qiime already exist; please, delete them first!"
  exit 1
fi

create_dir $INSTALL_DIR
download_archive
build

chmod -R ug+rwX,o+rX-w $INSTALL_DIR/$SOFTWARE_DIR

store_archive $ARCHIVE

# Deploy module
create_dir $MODULE_DIR
# Surround variable with "" since it contains a multiline text
module_file > $MODULE_DIR/${VERSION}_qiime
# Default module version file
echo "\
#%Module1.0
set ModulesVersion \"${VERSION}_qiime\"" > $MODULE_DIR/.version

chmod ug+rwX,o+rX-w $MODULE_DIR/${VERSION}_qiime $MODULE_DIR/.version

# Clean up temporary installation files if any
cd
rm -rf $INSTALL_DOWNLOAD

echo
echo "$SOFTWARE version ${VERSION}_qiime has been successfully installed in \$$INSTALL_HOME"
if [[ ! $INSTALL_HOME == 'MUGQIC_INSTALL_HOME' ]]
then
  echo "To install module in production, type '$0 MUGQIC_INSTALL_HOME' (no '\$' before parameter)"
fi

################################################################

sleep 30

################################################################
# Now install the all python libraries & dependencies for qiime
# Call the module to installall thel python libraries & dependencies for qiime
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/python_lib.2.7.10_qiime.sh 
