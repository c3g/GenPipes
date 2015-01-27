#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

SOFTWARE="tablemaker"
VERSION="2.1.1"  ## TO BE MODIFIED WITH e.g. 2.2.28+, 3.0, 0.1.19, etc.
ARCHIVE="$SOFTWARE"_"$VERSION".Linux_x86_64.tar.gz  ## TO BE MODIFIED WITH SPECIFIC ARCHIVE
ARCHIVE_URL="http://files.figshare.com/1529779/tablemaker_"$VERSION".Linux_x86_64.tar.gz"  ## TO BE MODIFIED WITH SPECIFIC URL
SOFTWARE_DIR="$SOFTWARE-$VERSION.Linux_x86_64"  ## TO BE MODIFIED WITH SPECIFIC SOFTWARE DIRECTORY IF NECESSARY

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE  ## TO BE MODIFIED WITH SPECIFIC COMMAND

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
prepend-path    PATH                \$root/ ;
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
