#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=Smoove
VERSION=0.2.3
ARCHIVE=${SOFTWARE,}-${VERSION}
ARCHIVE_URL=https://github.com/brentp/${SOFTWARE,}/releases/download/v${VERSION}/${SOFTWARE,}
SOFTWARE_DIR=${SOFTWARE,}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD

  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
  cp $ARCHIVE $INSTALL_DIR/$SOFTWARE_DIR/${SOFTWARE,}

  chmod 775 $INSTALL_DIR/$SOFTWARE_DIR/${SOFTWARE,}
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
prereq mugqic/LUMPY-SV/master_20190208 mugqic/samtools/1.9 mugqic/gsort/0.0.6 mugqic/htslib/1.9 mugqic/SVTyper/0.7.0 mugqic/python/2.7.14 mugqic/mosdepth/0.2.4 mugqic/bcftools/1.9 mugqic/duphold/0.1.1
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
