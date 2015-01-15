#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=amos
VERSION=3.1.0
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=http://downloads.sourceforge.net/project/$SOFTWARE/$SOFTWARE/$VERSION/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  # AMOS hard-codes dependency paths
  module load mugqic/perl/5.18.2 mugqic/python/2.7.8 mugqic/MUMmer/3.23 mugqic/ucsc/20140212
  ./configure \
    PERL=`which perl` \
    PYTHON=`which python` \
    NUCMER=`which nucmer` \
    DELTAFILTER=`which delta-filter` \
    SHOWCOORDS=`which show-coords` \
    BLAT=`which blat` \
    --prefix=$INSTALL_DIR/$SOFTWARE_DIR \
    --with-qmake-qt4=/usr/lib64/qt4/bin/qmake
  make
  make install

  # DELTAFILTER and SHOWCOORDS are ignored by ./configure for minimus2 (bug?)
  perl -pi -e s,"DELTAFILTER\s+=.*,DELTAFILTER=`which delta-filter`,g" $INSTALL_DIR/$SOFTWARE_DIR/bin/minimus2
  perl -pi -e s,"SHOWCOORDS\s+=.*,SHOWCOORDS=`which show-coords`,g" $INSTALL_DIR/$SOFTWARE_DIR/bin/minimus2
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
# IMPORTANT: to install amos in production, remove '../' before install_module.sh and move this script in parent dir
source $MODULE_INSTALL_SCRIPT_DIR/../install_module.sh $@
