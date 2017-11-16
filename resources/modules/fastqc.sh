#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=fastqc
VERSION=0.11.6.devel
#VERSION=0.11.5
ARCHIVE=$SOFTWARE-$VERSION.zip
ARCHIVE_URL=http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${VERSION}.zip
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE
  mv FastQC $SOFTWARE_DIR
  chmod +x $SOFTWARE_DIR/fastqc
  chmod -R g+w $SOFTWARE_DIR
  sed -i s,"#\!/usr/bin/perl.*,#\!/usr/bin/env perl,g" $SOFTWARE_DIR/fastqc
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

