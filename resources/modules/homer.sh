#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#
# HOMER
#

# Install the basic HOMER software. HOMER will be installed in the same directory that you place the configureHomer.pl program.  
# configureHomer.pl will attempt to check for required utilities and alert you to missing programs.
# To install packages (Genomes), simply use the –install option and the name(s) of the package(s).
# perl configureHomer.pl –install mouse (to download the mouse promoter set)
# perl configureHomer.pl –install mm8    (to download the mm8 version of the mouse genome)
# perl configureHomer.pl –install hg19r    (to download the hg19 repeat masked version of the human genome)

SOFTWARE=homer
VERSION=4.10
ARCHIVE_URL=http://homer.salk.edu/$SOFTWARE/configureHomer.pl
ARCHIVE=configureHomer.pl
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
  cp configureHomer.pl $INSTALL_DIR/$SOFTWARE_DIR/
  cd $INSTALL_DIR/$SOFTWARE_DIR

  module load mugqic/perl/5.22.1
  #module load mugqic/weblogo/3.3
  module load mugqic/ucsc/v346
  perl configureHomer.pl -install
  perl configureHomer.pl -install hg19
  perl configureHomer.pl -install hg38
  perl configureHomer.pl -install mm10
  perl configureHomer.pl -install mm9
  perl configureHomer.pl -install rn5
  perl configureHomer.pl -install rn6
  perl configureHomer.pl -install dm3

  # Update Perl scripts shebang
  find . -name "*.pl" | while read f ; do sed -i s,"#\!/usr/bin/perl -w,#\!/usr/bin/env perl\nuse warnings;,g" $f ; done
  find . -name "*.pl" | while read f ; do sed -i s,"#\!/usr/bin/perl,#\!/usr/bin/env perl,g" $f ; done
}

module_file() {
# Module file
echo "#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
setenv          HOMER_HOME          \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
