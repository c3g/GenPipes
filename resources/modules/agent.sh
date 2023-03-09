#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

HOST=`hostname`;
DNSDOMAIN=`dnsdomainname`;

SOFTWARE=AGeNT
VERSION=3.0.6
ARCHIVE=${SOFTWARE}_${VERSION}.zip
# Fill the form and downloa the zip archive from 
# https://www.agilent.com/en/product/next-generation-sequencing/hybridization-based-next-generation-sequencing-ngs/ngs-software/agent-232879
# Then put it in the archive folder e.g. $MUGQIC_INSTALL_HOME/archive
ARCHIVE_URL=
SOFTWARE_DIR=$SOFTWARE-$VERSION 

build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE
  mv -i ${SOFTWARE,,}  $INSTALL_DIR/$SOFTWARE_DIR
  mv *.html $INSTALL_DIR/$SOFTWARE_DIR/
  mv *.md $INSTALL_DIR/$SOFTWARE_DIR/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          AGENT_HOME          \$root;
setenv          TRIMMER             \$root/lib/trimmer-3.0.5.jar;
setenv          CREAKS              \$root/lib/creaks-1.0.5.jar;
setenv          LOCATIT             \$root/lib/locatit-2.0.5.jar
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
