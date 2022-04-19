#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

export NOPATCH=1
export NOWRAP=1

SOFTWARE=pandoc
VERSION=2.16.2
ARCHIVE=${SOFTWARE}-${VERSION}-linux-amd64.tar.gz
ARCHIVE_URL=https://github.com/jgm/pandoc/releases/download/${VERSION}/${ARCHIVE}
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar -xvf $ARCHIVE

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
prepend-path    PATH                \$root/bin
prepend-path    MANPATH             \$root/share/man/man1/
setenv          LANG                en_US.utf8
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
