#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=nextflow
VERSION=22.10.6
ARCHIVE=$SOFTWARE-$VERSION
ARCHIVE_URL=https://github.com/nextflow-io/${SOFTWARE}/releases/download/v${VERSION}/${SOFTWARE}-${VERSION}-all
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR/
  cp $ARCHIVE $INSTALL_DIR/$SOFTWARE_DIR/$SOFTWARE
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
setenv          NXF_OPTS            -Djdk.lang.Process.launchMechanism=vfork
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
