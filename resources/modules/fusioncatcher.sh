#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=FusionCatcher
VERSION=1.33
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/ndaniel/${SOFTWARE,,}/archive/refs/tags/${VERSION}.tar.gz
# ARCHIVE_URL=https://sourceforge.net/projects/${SOFTWARE,,}/files/${SOFTWARE,,}_v${VERSION}.zip
SOFTWARE_DIR=$SOFTWARE-$VERSION
PYTHON_VERSION=2.7.18
PYTHON_SHORT_VERSION=${PYTHON_VERSION%.*}

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  module load mugqic/python/$PYTHON_VERSION
  cd ${SOFTWARE,,}-$VERSION/tools
  ./install_tools.sh

  ## Do Not build data in the software folder
  #cd ../data
  #./download-human-db.sh

  # example of how to build the data
  # fusioncatcher-build  -g  homo_sapiens  -o /cvmfs/soft.mugqic/CentOS6/genomes/FusionCatcher/human_v108 

  cd $INSTALL_DOWNLOAD
  mv ${SOFTWARE,,}-$VERSION $INSTALL_DIR/$SOFTWARE_DIR
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
