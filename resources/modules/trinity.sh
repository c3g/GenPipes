#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=trinity
VERSION=2.2.0
ARCHIVE=${SOFTWARE}rnaseq-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/trinityrnaseq/trinityrnaseq/archive/v$VERSION.tar.gz
SOFTWARE_DIR=${SOFTWARE}rnaseq_$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  make -j12
  if [[ ${VERSION:0:1} == 2 ]]; then echo $VERSION; make -j12 plugins; fi

  # Install software
  cd $INSTALL_DOWNLOAD

  cd $SOFTWARE_DIR
  find . -name "*.pl" | while read f ; do sed -i s,"#\!/usr/bin/perl -w,#\!/usr/bin/env perl\nuse warnings;,g" $f ; sed -i s,"#\!/usr/bin/perl,#\!/usr/bin/env perl\nuse warnings;,g" $f ; done
  cd ..

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
