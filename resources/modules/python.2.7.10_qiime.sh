#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=python
VERSION=2.7.10_qiime
SETUPTOOLS_VERSION=20.9.0
# Remove the version last number
LIBVERSION=2.7
# Uppercase first P in python
ARCHIVE=${SOFTWARE^}-$VERSION.tgz
ARCHIVE_URL=http://www.python.org/ftp/$SOFTWARE/2.7.10/${SOFTWARE^}-2.7.10.tgz
SOFTWARE_DIR=${SOFTWARE^}-2.7.10

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  mv $SOFTWARE_DIR ${SOFTWARE^}-$VERSION
  SOFTWARE_DIR=${SOFTWARE^}-$VERSION

  cd $SOFTWARE_DIR
  # Compile with --enable-unicode=ucs4 to fix error "ImportError: numpy-1.8.1-py2.7-linux-x86_64.egg/numpy/core/multiarray.so: undefined symbol: PyUnicodeUCS2_AsASCIIString"
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR --enable-unicode=ucs4
  make -j12
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
setenv          PYTHONHOME          \$root
prepend-path    PATH                \$root/bin
prepend-path    MANPATH             \$root/share/man
prepend-path    LIBRARY_PATH        \$root/lib/
prepend-path    LD_LIBRARY_PATH     \$root/lib/
prepend-path    CPATH               \$root/include:\$root/include/python$LIBVERSION
prepend-path    PYTHONPATH          \$root/lib/python$LIBVERSION/site-packages:\$root/lib/python$LIBVERSION
setenv          QIIME_HOME          \$root/bin
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
else
  INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV
  MODULE_HOME=muqgic_dev
fi

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE
ARCHIVE_DIR=${!INSTALL_HOME}/archive

# Set module directory path by lowercasing $INSTALL_HOME and removing '_install_home' in it
MODULE_DIR=${!INSTALL_HOME}/modulefiles/`echo ${INSTALL_HOME,,} | sed 's/_install_home//'`/$SOFTWARE

echo "Installing $SOFTWARE version ${VERSION} in \$$INSTALL_HOME..."
echo

# Abort if software and/or module are already installed
if [[ -e $INSTALL_DIR/${SOFTWARE_DIR}_qiime || -e $MODULE_DIR/$VERSION ]]
then
  echo "$INSTALL_DIR/$SOFTWARE_DIR_qiime and/or $MODULE_DIR/$VERSION already exist; please, delete them first!"
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
module_file > $MODULE_DIR/${VERSION}
# Default module version file
echo "\
#%Module1.0
set ModulesVersion \"${VERSION}\"" > $MODULE_DIR/.version

chmod ug+rwX,o+rX-w $MODULE_DIR/${VERSION} $MODULE_DIR/.version

# Clean up temporary installation files if any
cd
rm -rf $INSTALL_DOWNLOAD

echo
echo "$SOFTWARE version ${VERSION} has been successfully installed in \$$INSTALL_HOME"
if [[ ! $INSTALL_HOME == 'MUGQIC_INSTALL_HOME' ]]
then
  echo "To install module in production, type '$0 MUGQIC_INSTALL_HOME' (no '\$' before parameter)"
fi

################################################################

sleep 30

################################################################

# Now installing Python libraries and qiime dependencies


#pip
echo "Installing pip..."
$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install pip
echo "pip installation complete !"

# cython
echo "Installing cython..."
$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install cython
echo "cython installation complete !"

# numpy
echo "Installing numpy..."
$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install numpy
python -c 'import numpy; print numpy.__version__, numpy.__file__'
echo "numpy installation complete !"

# biopython
echo "Installing biopython..."
$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install biopython
python -c 'import Bio; print Bio.__version__, Bio.__file__'
echo "biopython installation complete !"

# matplotlib requires dateutil and pyparsing dependencies
echo "Installing python-dateutil..."
$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install python-dateutil
echo "python-dateutil installation complete !"

echo "Installing pyparsing..."
$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install pyparsing
echo "pyparsing installation complete !"

echo "Installing matplotlib v1.4.3..."
$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install matplotlib==1.4.3
python -c 'import matplotlib; print matplotlib.__version__, matplotlib.__file__'
echo "matplotlib v1.4.3 installation complete !"

# HTseq
echo "Installing HTseq..."
$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install HTSeq
python -c 'import HTSeq; print HTSeq.__version__, HTSeq.__file__'
echo "HTseq installation complete !"

# bedtools-python has no version (version from doc: 0.1.0): install from master
echo "Installing bedtools-python..."
$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install https://github.com/arq5x/bedtools-python/archive/master.zip
python -c 'import bedtools; print  bedtools.__file__'
echo "bedtools-python installation complete !"

# PyVCF
echo "Installing PyVCF..."
$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install PyVCF
python -c 'import vcf; print vcf.__file__'
echo "PyVCF installation complete !"

# pysam
echo "Installing pysam..."
$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install pysam
python -c 'import pysam; print pysam.__version__, pysam.__file__'
echo "pysam installation complete !"

# nextworkx
echo "Insalling nextworkx..."
$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install networkx
python -c 'import networkx; print networkx.__version__, networkx.__file__'
echo "nextworkx installation complete !"

#qiime
echo "Insalling Qiime..."
$INSTALL_DIR/$SOFTWARE_DIR/bin/pip install qiime
echo "Qiime installation complete !"

# Add permissions
echo "adding the right permissions..."
chmod -R ug+rwX,o+rX-w $INSTALL_DIR/$SOFTWARE_DIR


echo "Python libraries and Qiime dependencies were successfully installed !!"

