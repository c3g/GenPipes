#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=perl
VERSION=5.28.2
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=http://www.cpan.org/src/5.0/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION
NOWRAP=1

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  ./Configure -des -Dusethreads -Dprefix=$INSTALL_DIR/$SOFTWARE_DIR
  make -j12
  make install
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE-$VERSION\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          PERL_HOME           \$root
prepend-path    PATH                \$root/bin
prepend-path    PERL5LIB            \$root/lib
prepend-path    PERL5LIB            \$root/lib/$VERSION
prepend-path    PERL5LIB            \$root/lib/site_perl/$VERSION
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
