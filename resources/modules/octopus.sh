#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=Octopus
VERSION=0.6.2-beta
ARCHIVE=${SOFTWARE,}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/luntergroup/${SOFTWARE,}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE,}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD

  # What follows is the installation from source (often problematic...)
  git clone --recursive https://github.com/luntergroup/${SOFTWARE,}.git -b v$VERSION

  mv ${SOFTWARE,} $SOFTWARE_DIR
  module load mugqic/python/3.6.5 mugqic/boost/1.66.0
  $SOFTWARE_DIR/scripts/install.py -cxx /cvmfs/soft.mugqic/yum/centos7/1.0/usr/local/bin/g++ -c /cvmfs/soft.mugqic/yum/centos7/1.0/usr/local/bin/gcc --boost $BOOST_ROOT --install-dependencies --download-forests

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
setenv          OCTOPUS_HOME        \$root
prepend-path    PATH                \$root/bin
prereq mugqic/boost/1.66.0
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
