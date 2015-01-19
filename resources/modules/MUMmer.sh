#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=MUMmer
VERSION=3.23
ARCHIVE=$SOFTWARE$VERSION.tar.gz
ARCHIVE_URL=http://downloads.sourceforge.net/project/mummer/mummer/$VERSION/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE --directory=$INSTALL_DIR/

  cd $INSTALL_DIR/$SOFTWARE_DIR
  make check
  make install

  # Update Perl script shebangs
  sed -i s,"#\!/usr/bin/perl.*,#\!/usr/bin/env perl,g" mapview mummerplot nucmer promer nucmer2xfig dnadiff
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
prepend-path    PERL5LIB            \$root/scripts/
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
