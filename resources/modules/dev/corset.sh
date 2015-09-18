#!/bin/bash
# Exit immediately on error
set -eu -o pipefail


SOFTWARE=Corset
VERSION=1.0.4  
ARCHIVE=version-1.04.zip
ARCHIVE_URL=https://github.com/Oshlack/Corset/archive/$ARCHIVE
SOFTWARE_DIR=Corset-version-1.04 

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD

build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE  ## TO BE MODIFIED WITH SPECIFIC COMMAND
  # Use samtools 
  mv -i $INSTALL_DOWNLOAD/$SOFTWARE_DIR $INSTALL_DIR/
  cd  $INSTALL_DIR/$SOFTWARE_DIR
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR --with-bam_inc=/cvmfs/soft.mugqic/CentOS6/software/samtools/samtools-1.0/ --with-bam_lib=/cvmfs/soft.mugqic/CentOS6/software/samtools/samtools-1.0/
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
prepend-path    PATH                \$root ;  
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
