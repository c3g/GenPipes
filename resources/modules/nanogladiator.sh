#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=NanoGLADIATOR
VERSION=1.0
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://sourceforge.net/projects/${SOFTWARE}/files/${SOFTWARE}_${VERSION}.tar.gz/download
SOFTWARE_DIR=${SOFTWARE}-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  mv ${SOFTWARE}_${VERSION} $SOFTWARE_DIR

  module load mugqic/R_Bioconductor/4.1.0_3.13
  cd $SOFTWARE_DIR/lib/F77
  R CMD SHLIB F4R.f
  R CMD SHLIB FastJointSLMLibrary.f

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
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
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
