#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=htslib
VERSION=1.12
ARCHIVE=$SOFTWARE-$VERSION.tar.bz2
ARCHIVE_URL=https://github.com/samtools/htslib/releases/download/${VERSION}/${ARCHIVE}
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar jxvf $ARCHIVE

  cd $SOFTWARE_DIR
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
prepend-path    LD_LIBRARY_PATH     \$root/lib
prepend-path    C_INCLUDE_PATH      \$root/include
prepend-path    CPP_INCLUDE_PATH    \$root/include
prepend-path    CPATH               \$root/include
prepend-path --delim \" \" LDFLAGS  -L\$root/lib
prepend-path --delim \" \" CPPFLAGS -I\$root/include
setenv          HTSLIB_HOME         \$root
setenv          HTSLIB_LIBRARY_DIR  \$root/lib
setenv          HTSLIB_INCLUDE_DIR  \$root/include
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
