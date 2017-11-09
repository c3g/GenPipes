#!/bin/bash
# Exit immediately on error
set -eu -o pipefail


SOFTWARE=Athena
VERSION=1.0
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/abishara/${SOFTWARE,}_meta/archive/${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE,}_meta-${VERSION}

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD

build() {
  cd $INSTALL_DOWNLOAD
  tar -zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  module load mugqic/python
  pip install .

  cd $INSTALL_DOWNLOAD
  mv $SOFTWARE_DIR $INSTALL_DIR/
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
module load mugqic/python mugqic/htslib mugqic/samtools mugqic/bwa mugqic/Canu mugqic/idba
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
