#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=CrossMap
VERSION=0.5.3
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
#ARCHIVE_URL=https://github.com/liguowang/CrossMap/archive/refs/tags/0.5.3.tar.gz
ARCHIVE_URL=https://sourceforge.net/projects/${SOFTWARE,,}/files/${VERSION}/version%20${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION
PYTHON_VERSION=3.9.1
PYTHON_SHORT_VERSION=${PYTHON_VERSION:0:3}
NOWRAP=1
NOPATCH=1

build() {
  cd $INSTALL_DOWNLOAD

  module load mugqic/python/$PYTHON_VERSION
#  pip install --prefix=$INSTALL_DIR/$SOFTWARE_DIR --ignore-installed ${SOFTWARE}==${VERSION}
  pip install --prefix=$INSTALL_DIR/$SOFTWARE_DIR --ignore-installed git+https://github.com/liguowang/CrossMap.git@${VERSION}

  # restting shebang so that CrossMap can be used with any version of Python
  for i in `find $INSTALL_DIR/$SOFTWARE_DIR/bin -type f`; do sed -i 's/^\#!.*/#!\/usr\/bin\/env python/' $i; done

  # Move the chain files (from the archive) to the CrossMap folder
  tar zxvf $ARCHIVE
  mv $SOFTWARE/data  $INSTALL_DIR/$SOFTWARE_DIR/

  # create a link of the python executable in the software bin folder
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
prepend-path    PYTHONPATH          \$root/lib/python$PYTHON_SHORT_VERSION
prepend-path    PYTHONPATH          \$root/lib/python${PYTHON_SHORT_VERSION}/site-packages
setenv          CROSSMAP_HOME       \$root
setenv          CROSSMAP_CHAINS     \$root/data
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

