#!/bin/bash
# Exit immediately on error
set -eu -o pipefail


SOFTWARE=AMP_Scanner
VERSION=V1_Proteome
ARCHIVE=${SOFTWARE}_${VERSION}.zip
ARCHIVE_URL=
# Use the zip file provided by Dan Veltri as there is no official release for now
# Refer to https://computationalgenomics.atlassian.net/browse/PRJBFX-1741
# File : https://drive.google.com/open?id=1DVFfOVJ_YANq5gxSOgakpjZD3xu9
SOFTWARE_DIR=${SOFTWARE}_${VERSION}


build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE

  sed -i -e 's/command -v bin/command -v $AMP_SCANNER_HOME\/bin/' $SOFTWARE_DIR/Proteome_AMP_Scanner.sh
  sed -i -e 's/ruby script/ruby $AMP_SCANNER_HOME\/script/' $SOFTWARE_DIR/Proteome_AMP_Scanner.sh
  sed -i -e 's/Rscript script/Rscript $AMP_SCANNER_HOME\/script/' $SOFTWARE_DIR/Proteome_AMP_Scanner.sh

  mv $SOFTWARE_DIR $INSTALL_DIR/
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
prepend-path    PATH                \$root/bin
setenv          AMP_SCANNER_HOME    \$root
prereq mugqic/Ruby/2.5.3 mugqic/R_Bioconductor/3.4.4_3.6 mugqic/emboss/6.6.0
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
