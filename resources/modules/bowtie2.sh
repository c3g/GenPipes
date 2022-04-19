#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bowtie2
VERSION=2.4.4
ARCHIVE=$SOFTWARE-$VERSION-source.zip
ARCHIVE_URL=https://sourceforge.net/projects/bowtie-bio/files/$SOFTWARE/$VERSION/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE

  cd $SOFTWARE_DIR
  make -j12
  make -j12 sra-deps
  make -j12 USE_SRA=1

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i ${SOFTWARE_DIR} $INSTALL_DIR/
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
