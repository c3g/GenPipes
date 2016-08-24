#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bowtie
VERSION=2.1.0
ARCHIVE=${SOFTWARE}2-$VERSION-source.zip
#ARCHIVE=$SOFTWARE-$VERSION-src.zip	# for version <2.0.0
#ARCHIVE_URL=https://sourceforge.net/projects/bowtie-bio/files/$SOFTWARE/$VERSION/$ARCHIVE	# for version <2.0.0
ARCHIVE_URL=https://sourceforge.net/projects/bowtie-bio/files/${SOFTWARE}2/$VERSION/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE}-$VERSION
#SOFTWARE_DIR=$SOFTWARE-$VERSION	# for version <2.0.0

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE

  mv ${SOFTWARE}2-$VERSION $SOFTWARE_DIR	# for version >= 2.0.0
  cd $SOFTWARE_DIR
  make

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
