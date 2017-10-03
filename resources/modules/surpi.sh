#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=SURPI
VERSION=1.0.18
ARCHIVE=${SOFTWARE,,}-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/chiulab/${SOFTWARE,,}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE,,}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE
  
  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR/

  cd $INSTALL_DIR/$SOFTWARE_DIR
  sed -i "s|#!/usr/bin/perl -w|#!/usr/bin/env perl\nuse warnings;|" *.pl
  sed -i "s|#!/usr/bin/perl|#!/usr/bin/env perl|" *.pl
  sed -i "s|#!/usr/bin/python|#!/usr/bin/env python|" *.py

  wget https://raw.github.com/attractivechaos/klib/master/khash.h
  wget http://chiulab.ucsf.edu/SURPI/software/fqextract.c
  gcc fqextract.c -o fqextract
  chmod +x fqextract

  gcc source/dropcache.c -o dropcache
  chmod u+s dropcache
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
module load mugqic/FastQValidator/0.1.1a
module load mugqic/AMOS/3.1.0
module load mugqic/ABySS/1.3.5
module load mugqic/RAPSearch2/2.12
module load mugqic/seqtk/1.0
module load mugqic/snap/0.15
module load mugqic/GenomeTools/1.5.9
module load mugqic/bio-playground/master
module load mugqic/python/2.7.13
module load mugqic/prinseq-lite/0.20.4
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

