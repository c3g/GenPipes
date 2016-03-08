#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

SOFTWARE=gemini 
VERSION=0.18.0  
ARCHIVE=${SOFTWARE}_v$VERSION.install.py
ARCHIVE_URL=https://raw.githubusercontent.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
#ARCHIVE_URL=https://raw.github.com/mbourgey/gemini/master/gemini/scripts/gemini_install.py
#ARCHIVE_URL=https://raw.github.com/mbourgey/$SOFTWARE/v$VERSION/$SOFTWARE/scripts/${SOFTWARE}_install.py 
SOFTWARE_DIR=$SOFTWARE-$VERSION 
PYTHON_VERSION=2.7.11

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  
  mkdir -p $SOFTWARE_DIR/shared_data $SOFTWARE_DIR/gemini_data
  module load mugqic/python/$PYTHON_VERSION
  python ${SOFTWARE}_v$VERSION.install.py $SOFTWARE_DIR  $SOFTWARE_DIR
  
  export PATH=$SOFTWARE_DIR/anaconda/bin:$PATH

  gemini update --dataonly --extra cadd_score --extra gerp_bp

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
setenv          GEMINI_BIN          \$root/bin
prepend-path    PATH                \$root/bin
prepend-path    PATH                \$root/anaconda/bin
prepend-path    PYTHONPATH          \$root/anaconda/lib/python2.7/site-packages
prepend-path    PYTHONPATH          \$root/anaconda/lib/python2.7
prepend-path    LD_LIBRARY_PATH     \$root/anaconda/lib/python2.7/site-packages
prepend-path    LD_LIBRARY_PATH     \$root/anaconda/lib/python2.7
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
