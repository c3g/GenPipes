#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################


# Java JDK (tested with OpenJDK 6) for GAGE
# Time::HiRes perl module for GeneMark-ES
# Boost (tested with v1.56.0) for Manta and E-MEM
# cmake (tested with v2.8.12) for Manta

module load cmake/3.3.1  jdk64/7u80 boost64/1.61.0  mugqic/perl/5.22.1 
#cpan --install Time::HighRes
SOFTWARE=quast  
VERSION=4.2
ARCHIVE="$SOFTWARE-$VERSION.tar.gz" 
ARCHIVE_URL="http://downloads.sourceforge.net/project/quast/$ARCHIVE" 
SOFTWARE_DIR=$SOFTWARE-$VERSION 
module load mugqic/python/2.7.8
# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE 

  cd $SOFTWARE_DIR
  python quast.py --test

  # Install software
  cd $INSTALL_DOWNLOAD  ## TO BE ADDED AND MODIFIED IF NECESSARY
  mv -i $SOFTWARE_DIR $INSTALL_DIR/  ## TO BE ADDED AND MODIFIED IF NECESSARY
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root ;  
prereq	mugqic/python
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
