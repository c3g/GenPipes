https://github.com/statgen/libStatGen/archive/v1.0.13.tar.gz
#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=libStatGen
VERSION=1.0.13
ARCHIVE=${SOFTWARE}-$VERSION.tar.gz 
ARCHIVE_URL=https://github.com/statgen/${SOFTWARE}/archive/v$VERSION.tar.gz
SOFTWARE_DIR=${SOFTWARE}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  cd $SOFTWARE_DIR
  make -j12

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
prepend-path    CPP_INCLUDE_PATH    \$root/include
prepend-path    LD_LIBRARY_PATH     \$root/
prepend-path    LIBRARY_PATH        \$root/
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
