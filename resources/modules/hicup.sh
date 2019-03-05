#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=HiCUP
VERSION=v0.7.2
ARCHIVE=${SOFTWARE,,}_${VERSION}.tar.gz
ARCHIVE_URL=https://www.bioinformatics.babraham.ac.uk/projects/${SOFTWARE,,}/${ARCHIVE}
SOFTWARE_DIR=${SOFTWARE,,}_$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  mv -i ${SOFTWARE}-master $INSTALL_DIR/$SOFTWARE_DIR

  cd $INSTALL_DIR/$SOFTWARE_DIR
  pwd
  ## change shebang to use loaded perl:
  sed -i "s|#!/usr/bin/perl -w|#!/usr/bin/env perl\nuse warnings;|" hicup*
  sed -i "s|#!/usr/bin/perl|#!/usr/bin/env perl|" hicup*

  if [ "`echo -e "${VERSION5/v/}\n${VERSION6/v/}" | sort -V  | head -n1`" = "0.5.9" ]
  then
    # Nothing to do if the version if 0.6.0 or above
    echo "version greater than 0.5.9" > /dev/null
  else
    # For early versions (earlier than 0.5.9), we change "-p1" for Bowtie2 to "-p8 --reorder" to force faster alignment in hicup_mapper in #Subroutine "map_file":
    echo "version smaller than 0.5.9" > /dev/null
    sed -i "s|\-\-very\-sensitive  \-x \$config{index} \-\-no\-unal \-p 1|-\-very\-sensitive  \-x \$config{index} \-\-no\-unal \-p 8 --reorder|" hicup_mapper 
  fi
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
