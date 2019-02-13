#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=MultiQC
VERSION=1.6
ARCHIVE=${SOFTWARE}-${VERSION}.zip
ARCHIVE_URL=https://github.com/ewels/$SOFTWARE/archive/v${VERSION}.zip
SOFTWARE_DIR=${SOFTWARE}-$VERSION
PYTHON_MODULE=mugqic/python/2.7.13

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {

  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE
#  mv ${SOFTWARE}-master $SOFTWARE_DIR

  mv -i $SOFTWARE_DIR $INSTALL_DIR/

  module load $PYTHON_MODULE

  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR/lib/python2.7/site-packages
  export PYTHONPATH=${PYTHONPATH}:$INSTALL_DIR/$SOFTWARE_DIR/lib/python2.7/site-packages
  cd $INSTALL_DIR/$SOFTWARE_DIR
  python setup.py install --prefix $INSTALL_DIR/$SOFTWARE_DIR
}


#Module definition to use
module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \"
       puts stderr \"\n\tVersion $VERSION\n\"
}

# only one version at a time
conflict multiqc

module-whatis \"MultiQC for for generating interactive analysis reports for pipeline\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
prepend-path    PYTHONPATH          \$root/lib/python2.7/site-packages
"
}


# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

