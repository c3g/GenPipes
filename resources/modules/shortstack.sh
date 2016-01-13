#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#
# Software_name ShortStack
#

SOFTWARE="ShortStack" 
#VERSION="2.1.0"
VERSION="3.3"
ARCHIVE=$VERSION.tar.gz
ARCHIVE_URL=http://github.com/MikeAxtell/ShortStack/archive/v$VERSION.tar.gz

SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # Install software
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
prepend-path    PATH                \$root
module load mugqic/bowtie/1.0.0 
module load mugqic_dev/ViennaRNA/1.8.3
module load mugqic/samtools/0.1.19
module load mugqic/ucsc/20140212
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
