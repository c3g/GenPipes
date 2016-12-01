#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=python
VERSION=2.7.12
SETUPTOOLS_VERSION=18.7.1
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

  EASY_INSTALL_PATH=$INSTALL_DIR/$SOFTWARE_DIR/bin/easy_install
  
  # pip
  ${EASY_INSTALL_PATH} pip
  PIP_PATH=$INSTALL_DIR/$SOFTWARE_DIR/bin/pip
  
  # cython
  #easy_install http://cython.org/release/Cython-0.23.4.tar.gz
  ${EASY_INSTALL_PATH} cython
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import cython; print cython.__version__, cython.__file__'

  # numpy
  ${EASY_INSTALL_PATH} numpy
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import numpy; print numpy.__version__, numpy.__file__'

  # biopython
  ${EASY_INSTALL_PATH} biopython
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import Bio; print Bio.__version__, Bio.__file__'

  # python-dateutil
  #easy_install http://labix.org/download/python-dateutil/python-dateutil-1.5.tar.gz
  ${EASY_INSTALL_PATH} python-dateutil
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import dateutil; print dateutil.__version__, dateutil.__file__'

  # pyparsing
  #easy_install https://sourceforge.net/projects/pyparsing/files/pyparsing/pyparsing-2.1.1/pyparsing-2.1.1.tar.gz/download
  ${EASY_INSTALL_PATH} pyparsing
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import pyparsing; print pyparsing.__version__, pyparsing.__file__'

  # matplotlib
  #easy_install https://github.com/matplotlib/matplotlib/archive/v1.5.1.tar.gz
  ${EASY_INSTALL_PATH} matplotlib
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import matplotlib; print matplotlib.__version__, matplotlib.__file__'

  # HTseq
  #easy_install https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.6.1p1.tar.gz
  ${EASY_INSTALL_PATH} HTseq
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import HTSeq; print HTSeq.__version__, HTSeq.__file__'

  # bedtools-python has no version (version from doc: 0.1.0): install from master
  ${EASY_INSTALL_PATH} https://github.com/arq5x/bedtools-python/archive/master.zip
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import bedtools; print  bedtools.__file__'

  # PyVCF
  #easy_install https://pypi.python.org/packages/source/P/PyVCF/PyVCF-0.6.8.tar.gz
  ${EASY_INSTALL_PATH} PyVCF
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import vcf; print vcf.__file__'

  # pysam
  #easy_install https://github.com/pysam-developers/pysam/archive/v0.9.0.tar.gz
  ${EASY_INSTALL_PATH} pysam
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import pysam; print pysam.__version__, pysam.__file__'

  # nextworkx
  ${EASY_INSTALL_PATH} https://github.com/networkx/networkx/archive/networkx-1.11.tar.gz
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import networkx; print networkx.__version__, networkx.__file__'

  # futures
  ${EASY_INSTALL_PATH} futures
  # future
  ${EASY_INSTALL_PATH} future
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import future; print future.__version__, future.__file__'

  # misopy
  ${EASY_INSTALL_PATH} misopy
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import misopy; print misopy.__version__, misopy.__file__'

  # qiime
  ${EASY_INSTALL_PATH} qiime
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import qiime; print qiime.__version__, qiime.__file__'

  # TEToolKit
  ${EASY_INSTALL_PATH} TEToolkit
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import TEToolkit; print TEToolkit.__file__'

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
