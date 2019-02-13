#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=samtools
VERSION=1.9
#VERSION=0.1.19
if [[ ${VERSION:0:1} == 1 ]]; then
  ARCHIVE=$SOFTWARE-$VERSION.tar.bz2
  ARCHIVE_URL=https://github.com/$SOFTWARE/$SOFTWARE/releases/download/$VERSION/$ARCHIVE
else
  ARCHIVE=$VERSION.tar.gz
  ARCHIVE_URL=https://github.com/$SOFTWARE/$SOFTWARE/archive/$ARCHIVE
fi
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  if [[ ${VERSION:0:1} == 1 ]]; then 
    tar jxvf $ARCHIVE
  else
    tar zxvf $ARCHIVE
  fi

  cd $SOFTWARE_DIR

  # Install software
  if [[ ${VERSION:0:1} == 1 ]]; then
    ./configure --enable-plugins --enable-libcurl prefix=$INSTALL_DIR/$SOFTWARE_DIR
    make -j12 all all-htslib
    make install install-htslib
    mv libbam.a $INSTALL_DIR/${SOFTWARE_DIR}/lib/
  else
    make -j12
    make -j12 razip
    mkdir -p bin
    cd bin
    ln -s ../misc/wgsim
    ln -s ../misc/md5sum-lite
    ln -s ../misc/md5fa
    ln -s ../misc/maq2sam-short
    ln -s ../misc/maq2sam-long
    for i in `ls ../misc/*.pl`; do ln -s $i; done
    for i in `ls ../misc/*.py`; do ln -s $i; done
    ln -s ../samtools
    ln -s ../razip
    cd ..
    cd $INSTALL_DOWNLOAD
    mv $SOFTWARE_DIR $INSTALL_DIR/

  fi

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
prepend-path    LIBRARY_PATH        \$root/lib
prepend-path    LD_LIBRARY_PATH     \$root/lib
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

