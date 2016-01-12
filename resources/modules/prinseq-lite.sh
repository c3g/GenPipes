#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=prinseq-lite
VERSION=0.20.4
#VERSION=0.20.3
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=http://sourceforge.net/projects/prinseq/files/standalone/$ARCHIVE/download
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  # Update Perl script shebangs
  sed -i s,"#\!/usr/bin/perl,#\!/usr/bin/env perl,g" prinseq-lite.pl prinseq-graphs.pl prinseq-graphs-noPCA.pl
  chmod +x prinseq-lite.pl prinseq-graphs.pl prinseq-graphs-noPCA.pl

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
