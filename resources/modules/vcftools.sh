#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=vcftools
<<<<<<< HEAD
VERSION=0.1.13
ARCHIVE=${SOFTWARE}_$VERSION.tar.gz
ARCHIVE_URL=http://sourceforge.net/projects/$SOFTWARE/files/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE}_$VERSION
=======
VERSION=0.1.14
ARCHIVE=${SOFTWARE}-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/$SOFTWARE/$SOFTWARE/releases/download/v${VERSION}/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE}-$VERSION
>>>>>>> e3c632527d6ee305d3d923bf8f6b50bb73329aca

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  FULL_PATH=$(readlink -f .)
  ./configure --prefix=$FULL_PATH
  make
  make install

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
prepend-path    PATH                \$root/bin
prepend-path    PERL5LIB            \$root/lib/perl5/site_perl
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
