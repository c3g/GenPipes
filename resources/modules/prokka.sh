#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE="Prokka"
VERSION="1.12"
ARCHIVE=${SOFTWARE,}-$VERSION.tar.gz
ARCHIVE_URL=http://www.vicbioinformatics.com/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE,}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR/

  MODULE_PERL=mugqic/perl/5.22.1
  MODULE_BLAST=mugqic/blast/2.3.0+
  cd $INSTALL_DIR/$SOFTWARE_DIR/bin
  module load $MODULE_PERL $MODULE_BLAST
  ./prokka --setupdb
  
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
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
