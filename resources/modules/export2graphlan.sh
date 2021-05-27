#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=export2graphlan
VERSION=0.22
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/SegataLab/${SOFTWARE,,}/archive/refs/tags/$VERSION.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION
PYTHON_VERSION=2.7.16
PYTHON_SHORT_VERSION=${PYTHON_VERSION:0:3}
NOWRAP=1
NOPATCH=1

build() {
  cd $INSTALL_DOWNLOAD
  module load mugqic/python/$PYTHON_VERSION
  pip install --prefix=$INSTALL_DIR/$SOFTWARE_DIR --ignore-installed ${SOFTWARE,,}==${VERSION}
  pip install --prefix=$INSTALL_DIR/$SOFTWARE_DIR hclust2
  ln -s $(which python) $INSTALL_DIR/$SOFTWARE_DIR/bin/python
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
prepend-path    PYTHONPATH          $PYTHONPATH
prepend-path    PYTHONPATH          \$root/lib/python${PYTHON_SHORT_VERSION}
prepend-path    PYTHONPATH          \$root/lib/python${PYTHON_SHORT_VERSION}/site-packages
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

