#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=circos-tools
VERSION=0.23
ARCHIVE=${SOFTWARE}-$VERSION.tgz
ARCHIVE_URL=http://circos.ca/distribution/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE}-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # Install software
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
prepend-path    PATH                \$root/tools

set circostools [glob -directory \$root/tools/ -- "*"]
foreach tool \$circostools {
  prepend-path  PATH                \$tool/bin
}
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
