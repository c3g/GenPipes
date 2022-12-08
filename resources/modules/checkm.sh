#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=CheckM
VERSION=1.2.1
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/Ecogenomics/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION
PYTHON_VERSION=3.10.4
PYTHON_SHORT_VERSION=${PYTHON_VERSION%.*}

build() {
  cd $INSTALL_DOWNLOAD
  module load mugqic/python/$PYTHON_VERSION
  pip3 install --prefix=$INSTALL_DIR/$SOFTWARE_DIR --ignore-installed numpy
  pip3 install --prefix=$INSTALL_DIR/$SOFTWARE_DIR --ignore-installed matplotlib
  pip3 install --prefix=$INSTALL_DIR/$SOFTWARE_DIR --ignore-installed pysam
  pip3 install --prefix=$INSTALL_DIR/$SOFTWARE_DIR --ignore-installed checkm-genome==$VERSION
  ln -s $(which python) $INSTALL_DIR/$SOFTWARE_DIR/bin/python
  ln -s $(which python3) $INSTALL_DIR/$SOFTWARE_DIR/bin/python3
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
setenv          CHECKM_DATA_PATH    /cvmfs/soft.mugqic/root/genomes/checkm_data
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

