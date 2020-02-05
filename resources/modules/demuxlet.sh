<<<<<<< HEAD
#!/bin/bash
=======
#!/usr/bin/env bash
>>>>>>> 1f66c6c99a2d466df7259210d21d3c202295b58f
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=demuxlet
VERSION=master_20190913
ARCHIVE=${SOFTWARE}-$VERSION.zip
ARCHIVE_URL=https://github.com/statgen/demuxlet/archive/master.zip
SOFTWARE_DIR=${SOFTWARE}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD

  # first install htslib in the same parent folder as demuxlet
  git clone https://github.com/samtools/htslib.git
  cd htslib
  autoheader     # If using configure, generate the header template...
  autoconf       # ...and configure script (or use autoreconf to do both)
  ./configure    # Optional but recommended, for choosing extra functionality
  make
  make install

  # Now install demuxlet
  cd ..  
  git clone --recursive https://github.com/statgen/demuxlet.git
  cd demuxlet
  autoreconf -vfi
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR
  make
  make install
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
