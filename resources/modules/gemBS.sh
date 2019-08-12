#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=gemBS 
VERSION=3.2.13
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/heathsc/${SOFTWARE}/archive/v3.2.0.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION
PYTHON_MODULE=mugqic/python/3.6.5

build() {
  cd $INSTALL_DOWNLOAD
  git clone --recursive https://github.com/heathsc/gemBS.git
  mv $SOFTWARE $SOFTWARE_DIR
  cd $SOFTWARE_DIR

  module load $PYTHON_MODULE
  python3 setup.py install --prefix $INSTALL_DIR/$SOFTWARE_DIR
  module unload $PYTHON_MODULE

  # resetting shebangs so that gemBS can be used with any version of Python
  sed -i 's/^\#!.*/#!\/usr\/bin\/env python/' $INSTALL_DIR/$SOFTWARE_DIR/bin/* 
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
prepend-path    PYTHONPATH          \$root/lib/python3.6/site-packages
puts stderr \"WARNING : the following modules needs python3 to run properly.\"
puts stderr \"Loading mugqic/python/3.6.5 would fit perfectly with this requirement.\"
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
