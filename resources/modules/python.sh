#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=python
VERSION=2.7.14
SETUPTOOLS_VERSION=36.4.0
# Remove the version last number
LIBVERSION=${VERSION%.[0-9]*}
# Uppercase first P in python
ARCHIVE=${SOFTWARE^}-$VERSION.tgz
ARCHIVE_URL=https://www.python.org/ftp/$SOFTWARE/$VERSION/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE^}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  # Compile with --enable-unicode=ucs4 to fix error "ImportError: numpy-1.8.1-py2.7-linux-x86_64.egg/numpy/core/multiarray.so: undefined symbol: PyUnicodeUCS2_AsASCIIString"
  #./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR --enable-unicode=ucs4 --with-zlib-dir=/usr/lib64 --with-ensurepip=install
  #./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR --enable-unicode=ucs4 --with-zlib-dir=$LFS/tools/usr/lib --with-ensurepip=install
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR --enable-unicode=ucs4 --with-zlib-dir=$LFS/tools/usr/lib --with-ensurepip=install
  make -j8
  make install

  echo "General Python installation done.... processing packages"

  # First update pip & setuptools
  PIP_PATH=$INSTALL_DIR/$SOFTWARE_DIR/bin/pip
  $PIP_PATH install --upgrade pip setuptools 

  # Then install the packages
  # cython
  #easy_install http://cython.org/release/Cython-0.23.4.tar.gz
  $PIP_PATH install --upgrade cython
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import cython; print cython.__version__, cython.__file__'

  # docopt
  $PIP_PATH install --upgrade docopt
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import docopt; print docopt.__version__, docopt.__file__'

  # numpy
  $PIP_PATH install --upgrade numpy
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import numpy; print numpy.__version__, numpy.__file__'

  # biopython
  $PIP_PATH install --upgrade biopython
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import Bio; print Bio.__version__, Bio.__file__'

  # python-dateutil
  #easy_install http://labix.org/download/python-dateutil/python-dateutil-1.5.tar.gz
  $PIP_PATH install --upgrade python-dateutil
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import dateutil; print dateutil.__version__, dateutil.__file__'

  # pyparsing
  #easy_install https://sourceforge.net/projects/pyparsing/files/pyparsing/pyparsing-2.1.1/pyparsing-2.1.1.tar.gz/download
  $PIP_PATH install --upgrade pyparsing
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import pyparsing; print pyparsing.__version__, pyparsing.__file__'

  # matplotlib
  #easy_install https://github.com/matplotlib/matplotlib/archive/v1.5.1.tar.gz
  $PIP_PATH install --upgrade matplotlib
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import matplotlib; print matplotlib.__version__, matplotlib.__file__'

  # HTseq
  #easy_install https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.6.1p1.tar.gz
  $PIP_PATH install --upgrade HTseq
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import HTSeq; print HTSeq.__version__, HTSeq.__file__'

  # bedtools-python has no version (version from doc: 0.1.0): install from master
  $PIP_PATH install --upgrade  https://github.com/arq5x/bedtools-python/archive/master.zip
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import bedtools; print bedtools.__file__'

  # PyVCF
  #easy_install https://pypi.python.org/packages/source/P/PyVCF/PyVCF-0.6.8.tar.gz
  $PIP_PATH install --upgrade PyVCF
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import vcf; print vcf.__file__'

  # pysam
  #easy_install https://github.com/pysam-developers/pysam/archive/v0.9.0.tar.gz
  $PIP_PATH install --upgrade pysam
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import pysam; print pysam.__version__, pysam.__file__'

  # networkx
  $PIP_PATH install --upgrade networkx
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import networkx; print networkx.__version__, networkx.__file__'

  # futures
  $PIP_PATH install --upgrade futures
  # future
  $PIP_PATH install --upgrade future
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import future; print future.__version__, future.__file__'

  # misopy
  $PIP_PATH install --upgrade misopy
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import misopy; print misopy.__version__, misopy.__file__'

  # qiime
  $PIP_PATH install --upgrade qiime
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import qiime; print qiime.__version__, qiime.__file__'

  # TEToolKit
  $PIP_PATH install --upgrade TEToolkit
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import TEToolkit; print TEToolkit.__file__'

  # deepTools
  $PIP_PATH install --upgrade deeptools
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import deeptools; print deeptools.__file__'

  # sklearn
  $PIP_PATH install --upgrade sklearn
  $INSTALL_DIR/$SOFTWARE_DIR/bin/python -c 'import sklearn; print sklearn.__version__; print sklearn.__file__'

  # RSeQC
  $PIP_PATH install --upgrade RSeQC
  
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

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

