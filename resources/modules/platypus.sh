#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=Platypus
VERSION=0.8.1

#echo "Prior to install the Patypus, you must download the archive $ARCHIVE manually, if not done already, from 
#Once downloaded, copy it in \$MUGQIC_INSTALL_HOME_DEV/archive/ or \$MUGQIC_INSTALL_HOME/archive/"

ARCHIVE=$SOFTWARE-$VERSION.tgz
ARCHIVE_URL=https://github.com/ekg/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE}_${VERSION}

htslib_module=mugqic/htslib/1.8

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  module load $htslib_module

  cd $SOFTWARE_DIR
  ./buildPlatypus.sh

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR ${INSTALL_DIR}/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
 puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
