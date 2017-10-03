#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=AffymetrixApt
VERSION=1.20.5
ARCHIVE=${SOFTWARE}-${VERSION}_x86_64_intel_linux.zip
# Archive has be manually downloaded form ThermoFisher website
# https://www.thermofisher.com/ca/en/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html
# THEN renamed as above and put in the $MUGQIC_INSTALL_HOME/archive or $MUGQIC_INSTALL_HOME_DEV/archive folder
ARCHIVE_URL=
SOFTWARE_DIR=${SOFTWARE}-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE
  mv apt-* $SOFTWARE_DIR
  
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
  chmod a+x $INSTALL_DIR/$SOFTWARE_DIR/bin/*
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"MUGQIC = Affymetrix $SOFTWARE toolkit\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
