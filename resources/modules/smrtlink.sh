#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=SMRTLink
VERSION=5.1.0
ARCHIVE=${SOFTWARE,,}-${VERSION}.26412.zip
ARCHIVE_URL=https://downloads.pacbcloud.com/public/software/installers/${ARCHIVE/-/_}
SOFTWARE_DIR=$SOFTWARE-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE
  ./smrtlink_5.1.0.26412.run --rootdir $INSTALL_DIR/$SOFTWARE_DIR --smrttools-only
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/smrtcmds/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

