#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=ivar
VERSION=1.2.2
ARCHIVE=$SOFTWARE-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/andersen-lab/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar -xzvf $ARCHIVE

  cd $SOFTWARE_DIR
  module load mugqic/htslib/1.10.2
  ./autogen.sh
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
puts stderr \"WARNING : iVar needs both samtools and htslib to be accessible in the envoronement.\"
puts stderr \"We recommend loading mugqic/samtools/1.10 and mugqic/htslib/1.10.2 \"
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
