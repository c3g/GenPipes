#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

SOFTWARE="sortmerna"
VERSION="2.0"
ARCHIVE="$VERSION.tar.gz"
ARCHIVE_URL="https://github.com/biocore/$SOFTWARE/archive/$VERSION.tar.gz"
SOFTWARE_DIR="$SOFTWARE-$VERSION"



build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE  ## TO BE MODIFIED WITH SPECIFIC COMMAND

  cd $SOFTWARE_DIR

  bash ./build.sh

  # Install software
  cd $INSTALL_DOWNLOAD  ## TO BE ADDED AND MODIFIED IF NECESSARY
  mv -i $SOFTWARE_DIR $INSTALL_DIR/  ## TO BE ADDED AND MODIFIED IF NECESSARY
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root;
setenv          SORTMERNA_HOME        \$root;
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
