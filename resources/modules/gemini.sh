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
VERSION=0.18.3  
ARCHIVE=${SOFTWARE}_v$VERSION.install.py
ARCHIVE_URL=https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
#ARCHIVE_URL=https://raw.githubusercontent.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
#ARCHIVE_URL=https://raw.github.com/mbourgey/gemini/master/gemini/scripts/gemini_install.py
#ARCHIVE_URL=https://raw.github.com/mbourgey/$SOFTWARE/v$VERSION/$SOFTWARE/scripts/${SOFTWARE}_install.py 
SOFTWARE_DIR=$SOFTWARE-$VERSION 
PYTHON_VERSION=2.7.11

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DIR
  mv tmp/* .

#  mkdir -p $SOFTWARE_DIR/shared_data $INSTALL_HOME/software/${SOFTWARE}/gemini_data
  mkdir -p $SOFTWARE_DIR/shared_data $INSTALL_DIR/shared
  module load mugqic/python/$PYTHON_VERSION
#  python ${SOFTWARE}_v$VERSION.install.py $SOFTWARE_DIR $INSTALL_HOME/software/${SOFTWARE}/gemini_data
  python ${SOFTWARE}_v$VERSION.install.py $SOFTWARE_DIR $INSTALL_DIR/shared

  export PATH=$INSTALL_DIR/shared/anaconda/bin:$PATH

  gemini update --dataonly --extra cadd_score --extra gerp_bp

  # Install software
#  cd $INSTALL_DOWNLOAD  ## TO BE ADDED AND MODIFIED IF NECESSARY
#  mv -i $SOFTWARE_DIR $INSTALL_DIR/  ## TO BE ADDED AND MODIFIED IF NECESSARY
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
set             anaconda_root       $INSTALL_DIR/shared/anaconda
setenv          GEMINI_BIN          \$root/bin
prepend-path    PATH                \$root/bin
prepend-path    PATH                \$anaconda_root/bin
prepend-path    PYTHONPATH          \$anaconda_root/lib/python2.7/site-packages
prepend-path    PYTHONPATH          \$anaconda_root/lib/python2.7
prepend-path    LD_LIBRARY_PATH     \$anaconda_root/lib/python2.7/site-packages
prepend-path    LD_LIBRARY_PATH     \$anaconda_root/lib/python2.7
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
