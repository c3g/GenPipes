#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bowtie2
VERSION=2.3.1
ARCHIVE=$SOFTWARE-$VERSION.zip
ARCHIVE_URL=https://sourceforge.net/projects/bowtie-bio/files/$SOFTWARE/$VERSION/${SOFTWARE}-${VERSION}-linux-x86_64.zip
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE

#  cd $SOFTWARE_DIR	# uncomment this if source code has been downloaded instead of binaries
#  make -j12

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
