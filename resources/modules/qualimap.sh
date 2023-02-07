#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

echo "Be sure to load java & R modules before using Qualimap"

SOFTWARE="qualimap"
VERSION="2.2.2dev"
ARCHIVE=$SOFTWARE-$VERSION.zip
#ARCHIVE_URL=https://bitbucket.org/kokonech/$SOFTWARE/downloads/${SOFTWARE}_v${VERSION}.zip
ARCHIVE_URL=https://bitbucket.org/kokonech/$SOFTWARE/downloads/${SOFTWARE}-build-11-11-19.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

R_MODULE=mugqic/R_Bioconductor/3.5.3_3.8

build() {
  cd $INSTALL_DOWNLOAD
#  unzip $ARCHIVE
  tar zxvf $ARCHIVE
#  mv ${SOFTWARE}_v${VERSION} $SOFTWARE_DIR
  mv ${SOFTWARE}-build-11-11-19 $SOFTWARE_DIR

  # Install software
  mv -i $SOFTWARE_DIR $INSTALL_DIR/

  # install required R-packages
  source /etc/profile.d/modules.sh
  module load $R_MODULE
  Rscript $INSTALL_DIR/$SOFTWARE_DIR/scripts/installDependencies.r
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
setenv          QUALIMAP_HOME       \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

