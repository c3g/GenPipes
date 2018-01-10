#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=gemini 
VERSION=0.20.1 
ARCHIVE=${SOFTWARE}_v$VERSION.install.py
ARCHIVE_URL=https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
SOFTWARE_DIR=$SOFTWARE-$VERSION 

build() {
  cd $INSTALL_DOWNLOAD

  # Install software but databases are not installed here.
  # Database installation has to be done manually afterwards to avoid any version conflict
  module load mugqic/python
  python $ARCHIVE --nodata $INSTALL_DIR/$SOFTWARE_DIR $INSTALL_DIR/$SOFTWARE_DIR/shared
  module unload mugqic/python

  echo ${VERSION:2:2}
  if [[ ! ${VERSION:2:2} > 18 ]]; then
    if [[ -f $ARCHIVE_DIR/$SOFTWARE-${VERSION}.tar.gz ]];
    then
      cp -a $ARCHIVE_DIR/$SOFTWARE-${VERSION}.tar.gz ./
      tar -zxvf $SOFTWARE-${VERSION}.tar.gz
    else
      wget --no-check-certificate https://github.com/arq5x/${SOFTWARE}/archive/v${VERSION}.tar.gz --output-document=$SOFTWARE-${VERSION}.tar.gz
      tar -zxvf $SOFTWARE-${VERSION}.tar.gz
      mv $SOFTWARE-${VERSION}.tar.gz $ARCHIVE_DIR/
    fi

    $INSTALL_DIR/$SOFTWARE_DIR/shared/anaconda/bin/pip install -r $SOFTWARE-${VERSION}/requirements.txt

    sed -i 's,$SOFTWARE-${VERSION}/,,' $INSTALL_DIR/$SOFTWARE_DIR/shared/gemini-config.yaml

  fi
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
set             anaconda_root       \$root/shared/anaconda
setenv          GEMINI_BIN          \$root/bin
prepend-path    PATH                \$root/bin
prepend-path    PATH                \$anaconda_root/bin
prepend-path    PYTHONHOME          \$anaconda_root
prepend-path    PYTHONPATH          \$anaconda_root/lib/python2.7/site-packages
prepend-path    PYTHONPATH          \$anaconda_root/lib/python2.7
prepend-path    LD_LIBRARY_PATH     \$anaconda_root/lib/python2.7/site-packages
prepend-path    LD_LIBRARY_PATH     \$anaconda_root/lib/python2.7
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
