#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

SOFTWARE=KronaTools
VERSION=2.6.1
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/marbl/${SOFTWARE:0:5}/releases/download/v${VERSION}/${SOFTWARE}-${VERSION}.tar
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar xvf $ARCHIVE

  # Move the forlder before installing to avoid creating wrong links
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
  cd $INSTALL_DIR/$SOFTWARE_DIR

  # Install software
  perl install.pl --prefix ./

  # Change the shebang of all the perl files in 'scripts/' directory
  cd scripts/
  for i in `ls *.pl`; do sed -ir 's/perl/env perl/' $i; done
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
prepend-path    PERL5LIB            \$root/lib
setenv          KRONA_HOME          \$root/scripts
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

