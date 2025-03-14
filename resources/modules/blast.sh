#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=blast
VERSION=2.11.0+
ARCHIVE=ncbi-$SOFTWARE-$VERSION-src.tar.gz
ARCHIVE_URL=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${VERSION%+}/$ARCHIVE
SOFTWARE_DIR=ncbi-$SOFTWARE-$VERSION


build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd ${SOFTWARE_DIR}-src/c++
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR
  make -j12
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
prepend-path    BLASTDB             \$::env(MUGQIC_INSTALL_HOME)/genomes/blast_db
setenv          BLAST_HOME          \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
