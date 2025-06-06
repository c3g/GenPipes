#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bcftools
VERSION=1.16
ARCHIVE=$SOFTWARE-$VERSION.tar.bz2
ARCHIVE_URL=https://github.com/samtools/bcftools/releases/download/${VERSION}/${ARCHIVE}
SOFTWARE_DIR=$SOFTWARE-$VERSION
NOWRAP=1

build() {
  cd $INSTALL_DOWNLOAD
  tar jxvf $ARCHIVE

  cd $SOFTWARE_DIR
  autoheader
  autoconf
  ./configure --enable-libgsl --enable-perl-filters
  make -j12
  # Install software
  make -j12 prefix=$INSTALL_DIR/${SOFTWARE_DIR} install
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
setenv          BCFTOOLS_PLUGINS    \$root/libexec/bcftools
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
