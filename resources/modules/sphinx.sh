#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=sphinx
VERSION=master
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
#ARCHIVE_URL=https://bitbucket.org/mugqic/$SOFTWARE/downloads/$ARCHIVE
ARCHIVE_URL=https://github.com/sphinx-doc/sphinx/archive/master.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # Install software
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
prepend-path    PATH                \$root/tools
prepend-path    PATH                \$root/perl-tools
prepend-path    PATH                \$root/R-tools
prepend-path    PATH                \$root/python-tools
prepend-path    PATH                \$root/RRNATagger-tools
prepend-path    PERL5LIB            \$root/perl-tools
setenv          R_TOOLS             \$root/R-tools
setenv          PERL_TOOLS          \$root/perl-tools
setenv          PYTHON_TOOLS        \$root/python-tools
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
