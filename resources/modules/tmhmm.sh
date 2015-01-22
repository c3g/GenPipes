#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

## NOTE (FL): 
# - Tools from cbs.dtu.dk under some sort of danish license and need to be downloaded manually.
# - Assumes gnuplot, gs and are avail on the system
# - Empty output problem: http://sourceforge.net/p/trinotate/mailman/message/31228028/
#
# You need an executable of the program decodeanhmm that runs under
# Unix. The program may already be in bin/decodeanhmm.
# 
# The scripts require perl 5.x
# 
# For plotting gnuplot is needed (making postscript plots).
# 
# When generating html output the postscript plots are converted to
# gif, and for this you need the programs ghostscript (gs) and ppmtogif.
# 
# After unpacking the directory you should
# 
# 1. Insert the correct path for perl 5.x in the first line of the scripts
#    bin/tmhmm and bin/tmhmmformat.pl (if not /usr/local/bin/perl).
# 2. Make sure you have an executable version of decodeanhmm in the bin
#    directory.
# 3. Include the directory containing tmhmm in your path.
# 4. Read the TMHMM2.0.guide.html.
# 5. Run the program.

SOFTWARE=tmhmm
VERSION=2.0c
ARCHIVE=$SOFTWARE-$VERSION.Linux.tar.gz
ARCHIVE_URL=http://www.dropbox.com/s/2lgbq3ci3zvvm2d/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR
  # ppm2gif, could not find this on the web. Link to ppm2tiff may work instead
  ln -s `which ppm2tiff` bin/ppm2gif

  # Update Perl script shebang
  sed -i "s,#!/usr/.*/perl.*,#!/usr/bin/env perl," bin/tmhmm bin/tmhmmformat.pl

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
prepend-path    PATH                \$root/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
