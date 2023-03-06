#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=vcf2db
VERSION=master.c391574_20210730
ARCHIVE=$SOFTWARE-${VERSION}.zip
ARCHIVE_URL=https://github.com/quinlan-lab/${SOFTWARE}/archive/master.zip
SOFTWARE_DIR=$SOFTWARE-${VERSION}
PYTHON_VERSION=3.6.5
PYTHON_SHORT_VERSION=${PYTHON_VERSION:0:3}
NOWRAP=1
NOPATCH=1

build() {
  cd $INSTALL_DIR
  git clone https://github.com/quinlan-lab/$SOFTWARE
  cd $SOFTWARE
  module load mugqic/python/$PYTHON_VERSION
  pip install --prefix=$INSTALL_DIR/$SOFTWARE_DIR --ignore-installed -r requirements.txt

  mv vcf2db.py $INSTALL_DIR/$SOFTWARE_DIR/bin/
  chmod a+rx $INSTALL_DIR/$SOFTWARE_DIR/bin/vcf2db.py

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
prepend-path    PYTHONPATH          \$root/lib/python${PYTHON_SHORT_VERSION}
prepend-path    PYTHONPATH          \$root/lib/python${PYTHON_SHORT_VERSION}/site-packages
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
