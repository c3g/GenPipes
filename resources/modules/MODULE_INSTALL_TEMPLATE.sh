#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

SOFTWARE=software_name  ## TO BE MODIFIED WITH e.g. blast, hmmer, samtools, etc.
VERSION=0.0.0  ## TO BE MODIFIED WITH e.g. 2.2.28+, 3.0, 0.1.19, etc.
ARCHIVE=$SOFTWARE-$VERSION.(tar.gz|zip|tar.bz2)  ## TO BE MODIFIED WITH SPECIFIC ARCHIVE
ARCHIVE_URL=http://www.software_lab.org/download/$ARCHIVE  ## TO BE MODIFIED WITH SPECIFIC URL
SOFTWARE_DIR=$SOFTWARE-$VERSION  ## TO BE MODIFIED WITH SPECIFIC SOFTWARE DIRECTORY IF NECESSARY

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  (tar zxvf|unzip|tar jxvf) $ARCHIVE  ## TO BE MODIFIED WITH SPECIFIC COMMAND

  cd $SOFTWARE_DIR
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR  ## TO BE ADDED AND MODIFIED IF NECESSARY
  make -j12  ## TO BE ADDED AND MODIFIED IF NECESSARY
  make install  ## TO BE ADDED IF NECESSARY

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
prepend-path    PATH                \$root/bin ;  ## TO BE ADDED IF NECESSARY
prepend-path    PATH                \$root/other_tools/bin ;  ## TO BE ADDED AND MODIFIED IF NECESSARY
setenv          ${SOFTWARE}_JAR     \$root/$SOFTWARE-$VERSION.jar ;  ## TO BE ADDED AND MODIFIED IF NECESSARY
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
