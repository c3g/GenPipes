#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

# WebLogo requires a recent version of ghostscript to create PNG and PDF output, and pdf2svg to generate SVG output.
# WebLogo version 3 is written in python. It is necessary to have python 2.5, 2.6 or 2.7 and the extension package numpy installed before WebLogo will run.

SOFTWARE=deepTools
VERSION=2.5.3
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/fidelram/${SOFTWARE}/archive/${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION
PYTHON_VERSION=2.7.13

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  cd $SOFTWARE_DIR
  module load mugqic/python/$PYTHON_VERSION
  python setup.py install --prefix $INSTALL_DIR/$SOFTWARE_DIR
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

