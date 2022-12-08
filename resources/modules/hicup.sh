#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=HiCUP
VERSION=0.8.3
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/StevenWingett/${SOFTWARE}/archive/refs/tags/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE}-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  mv -i $SOFTWARE_DIR $INSTALL_DIR/

  cd $INSTALL_DIR/$SOFTWARE_DIR

  ## change shebang to use loaded perl:
  sed -i "s|#!/usr/bin/perl -w|#!/usr/bin/env perl\nuse warnings;|" hicup*
  sed -i "s|#!/usr/bin/perl|#!/usr/bin/env perl|" hicup*

  # For early versions (earlier than 0.5.9), we change "-p1" for Bowtie2 to "-p8 --reorder" to force faster alignment in hicup_mapper in #Subroutine "map_file":
  # echo "version smaller than 0.5.9" > /dev/null
  # sed -i "s|\-\-very\-sensitive  \-x \$config{index} \-\-no\-unal \-p 1|-\-very\-sensitive  \-x \$config{index} \-\-no\-unal \-p 8 --reorder|" hicup_mapper
}

#Module definition to use
module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"HiCUP aligner for Hi-C\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
