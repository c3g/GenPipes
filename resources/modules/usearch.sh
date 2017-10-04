#!/bin/sh
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=usearch
VERSION=10.0.240
ARCHIVE=${SOFTWARE}${VERSION}
ARCHIVE_URL=
echo "Prior to install the usearch module, you must download the archive manually, if not done already, from http://www.drive5.com/usearch/download.html since it requires a license agreement.
Once downloaded, copy it in \$MUGQIC_INSTALL_HOME_DEV/archive/ or \$MUGQIC_INSTALL_HOME/archive/ and rename it as ${SOFTWARE}${VERSION}"
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  mkdir $INSTALL_DIR/$SOFTWARE_DIR
  echo "mv $INSTALL_DOWNLOAD/${SOFTWARE}* $INSTALL_DIR/$SOFTWARE_DIR/${SOFTWARE}${VERSION}"
  mv $INSTALL_DOWNLOAD/${SOFTWARE}* $INSTALL_DIR/$SOFTWARE_DIR/${SOFTWARE}${VERSION}
  cd $INSTALL_DIR/$SOFTWARE_DIR/ 
  ln -s ${SOFTWARE}${VERSION} ${SOFTWARE}
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE-$VERSION \" ;
}
module-whatis \"$SOFTWARE  \" ; 
                      
set             root                $INSTALL_DIR/$SOFTWARE_DIR ;
prepend-path    PATH                \$root ;  
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
