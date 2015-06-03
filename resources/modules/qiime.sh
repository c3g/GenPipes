#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

SOFTWARE=Qiime 
VERSION=1.9.1
ARCHIVE=$SOFTWARE-$VERSION.tar.gz  
ARCHIVE_URL=http://qiime.org/install/install.html
SOFTWARE_DIR=$SOFTWARE-$VERSION  

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  
  cd $INSTALL_DIR
  mkdir -p $SOFTWARE_DIR/lib/python2.7/site-packages
  # YOU MUST EXPORT the lib/python2.7/site-packages/ directory in your PYTHONPATH.
  export PYTHONPATH=$INSTALL_DIR/$SOFTWARE_DIR/lib/python2.7/site-packages:$PYTHONPATH
  easy_install --prefix=$INSTALL_DIR/$SOFTWARE_DIR/ qiime
  
  # Testing the QIIME installation
  $INSTALL_DIR/$SOFTWARE_DIR/bin/print_qiime_config.py -t
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR/bin
setenv          QIIME_HOME         \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
