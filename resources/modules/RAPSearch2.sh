#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=RAPSearch2
VERSION=2.12
ARCHIVE=${SOFTWARE%?}-$VERSION.tar.gz
ARCHIVE_URL=https://sourceforge.net/projects/${SOFTWARE,,}/files/`echo "${SOFTWARE%?}" | tr '[:upper:]' '[:lower:]'`${VERSION}_pair_64bits.tar.gz	# for version < 2.15
#ARCHIVE_URL=https://sourceforge.net/projects/${SOFTWARE,,}/files/${SOFTWARE%?}${VERSION}_64bits.tar.gz							# for version >= 2.15
SOFTWARE_DIR=${SOFTWARE}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  mv ${SOFTWARE%?}${VERSION}_pair_64bits $SOFTWARE_DIR		# for version < 2.12
#  mv ${SOFTWARE%?}${VERSION}_64bits $SOFTWARE_DIR		# for version >= 2.12

  cd $SOFTWARE_DIR
  ./install

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
prepend-path    PATH                \$root/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

