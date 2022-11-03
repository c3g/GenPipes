#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

REPO="rpackages" 
SOFTWARE=mugqic_R_packages
VERSION=1.0.7
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://bitbucket.org/mugqic/rpackages/downloads/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION
PCKGS="gqUtils gqSeqUtils"
R_MODULE="mugqic/R_Bioconductor"

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # Install software
  cd $SOFTWARE_DIR
  module load $R_MODULE
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
  R CMD INSTALL -l $INSTALL_DIR/$SOFTWARE_DIR $PCKGS 
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    RLIBS               \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
