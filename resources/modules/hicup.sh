#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=hicup
VERSION=v0.5.9
ARCHIVE=${SOFTWARE}_${VERSION}.tar.gz
ARCHIVE_URL=https://www.bioinformatics.babraham.ac.uk/projects/${SOFTWARE}/${ARCHIVE}
SOFTWARE_DIR=${SOFTWARE}_$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE


  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR/


  ## change shebang to use loaded perl:
  # sed -i "s|#!/usr/bin/perl -w|#!/usr/bin/env perl\nuse warnings;|" *.pl
  # sed -i "s|#!/usr/bin/perl|#!/usr/bin/env perl|" *.pl
  # sed -i "s|#!/usr/bin/python|#!/usr/bin/env python|" *.py

  cd $INSTALL_DIR/$SOFTWARE_DIR

  sed -i "s|#!/usr/bin/perl -w|#!/usr/bin/env perl\nuse warnings;|" hicup_*
  sed -i "s|#!/usr/bin/perl|#!/usr/bin/env perl|" hicup_*


  ## change "-p1" for Bowtie2 to "-p8 --reorder" to force faster alignment in hicup_mapper in #Subroutine "map_file":
  
  sed -i "s|\-\-very\-sensitive  \-x \$config{index} \-\-no\-unal \-p 1|-\-very\-sensitive  \-x \$config{index} \-\-no\-unal \-p 8 --reorder|" hicup_mapper 
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

