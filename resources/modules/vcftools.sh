#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=vcftools
VERSION=0.1.14
#VERSION=0.1.13
#VERSION=0.1.11
ARCHIVE=${SOFTWARE}-$VERSION.tar.gz
#ARCHIVE_URL=https://sourceforge.net/projects/$SOFTWARE/files/${SOFTWARE}_${VERSION}.tar.gz/download		# for version earlier than 0.1.14
ARCHIVE_URL=https://github.com/$SOFTWARE/$SOFTWARE/releases/download/v${VERSION}/$ARCHIVE			# for version 0.1.14
SOFTWARE_DIR=${SOFTWARE}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  FULL_PATH=$(readlink -f .)
  ./configure --prefix=$FULL_PATH			#  # for version 0.1.14
  make -j12
  make install

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
}

#Module definition to use for version 0.1.13
#module_file() {
#echo "\
##%Module1.0
#proc ModulesHelp { } {
#  puts stderr \"\tMUGQIC - $SOFTWARE \"
#}
#module-whatis \"$SOFTWARE\"
#
#set             root                $INSTALL_DIR/$SOFTWARE_DIR
#prepend-path    PATH                \$root/bin
#prepend-path    PERL5LIB            \$root/lib/perl5/site_perl
#"
#}

#Module definition to use for version 0.1.14
module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
prepend-path    PERL5LIB            \$root/share/perl5
"
}


# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
