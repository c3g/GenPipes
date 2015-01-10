#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=trinity
VERSION=20140413p1
ARCHIVE=${SOFTWARE}rnaseq_r$VERSION.tar.gz
ARCHIVE_URL=http://downloads.sourceforge.net/project/trinityrnaseq/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE}rnaseq_r$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  sed -i 's/cd trinity-plugins\/jellyfish \&\& \.\/configure.*/cd trinity-plugins\/jellyfish \&\& .\/configure CC=gcc CXX=g++ --enable-static --disable-shared --prefix=`pwd` \&\& $(MAKE) LDFLAGS=-all-static AM_CPPFLAGS="-Wall -Wnon-virtual-dtor -I"`pwd`/' Makefile
  # Disable "-" to "_" substitution in abundance matrix since we don't want sample names to be modified
  sed -ri 's/(\$file =\~ s\/\\-\/_\/g;)/#\1/' util/abundance_estimates_to_matrix.pl
  make

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
prepend-path    PATH                \$root
prepend-path    PATH                \$root/util
prepend-path    PATH                \$root/util/RSEM_util
prepend-path    PATH                \$root/Analysis/DifferentialExpression
setenv          TRINITY_HOME        \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
