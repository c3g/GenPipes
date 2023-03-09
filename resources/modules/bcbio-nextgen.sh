#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bcbio-nextgen
VERSION=1.2.9
ARCHIVE=$SOFTWARE-$VERSION.tar.gz 
ARCHIVE_URL=https://github.com/bcbio/${SOFTWARE}/archive/refs/tags/v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION 
PYTHON_VERSION=3.10.4
PYTHON_SHORT_VERSION=${PYTHON_VERSION%.*}

build() {
  module load mugqic/python/$PYTHON_VERSION
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  python3 $SOFTWARE_DIR/scripts/bcbio_nextgen_install.py $INSTALL_DIR/$SOFTWARE_DIR --tooldir=$INSTALL_DIR/$SOFTWARE_DIR/tools --nodata

}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
prepend-path    PATH                \$root/bin
prepend-path    PYTHONPATH          $PYTHONPATH
prepend-path    PYTHONPATH          \$root/lib/python${PYTHON_SHORT_VERSION}
prepend-path    PYTHONPATH          \$root/lib/python${PYTHON_SHORT_VERSION}/site-packages
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
