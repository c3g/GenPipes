#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

## NOTE (FL): Tools from cbs.dtu.dk under some sort of danish license and need to be downloaded manually.

# INSTALLATION
# 
# 1. Decide where in your file system you wish to keep the software.  Uncompress
#    and untar the package in that location:
# 
#         cat signalp-4.1.tar.Z | uncompress | tar xvf -
# 
#    This will produce a directory  'signalp-4.1'.  The size of the uncompressed
#    package will not exceed 30 Mb.
# 
# 2. Edit the paragraph labeled  "GENERAL SETTINGS, CUSTOMIZE ..." in the top of
#    the file 'signalp'. The following twovmandatory variables need to be set:
#    
#         SIGNALP         full path to the signalp-4.1 directory on your system
#         outputDir       where to store temporary files (writable to all users)
# 
#    In addition,  for practical reasons,  it is possible to limit the number of
#    input sequences allowed per run (MAX_ALLOWED_ENTRIES).
#    3. Test SignalP on the 10 eukaryotic sequences shipped with the package:
# 
#            > cd $SIGNALP
#            > ./signalp -t euk -f short test/euk10.fsa > euk10.fsa.short_out
#            > ./signalp -t euk -f longtest/euk10.fsa > euk10.fsa.long_out
#            > ./signalp -t euk -f all test/euk10.fsa > euk10.fsa.all_out
#            > ./signalp -t euk -f summary test/euk10.fsa > euk10.fsa.summary_out
# 
#       The output files "*_out"  should be identical to the corresponding files in
#       'signalp-4.1/test'.
# 
#    4. Move or copy the 'signalp' script to a directory in the users' path.
# 
#    5. Move or copy the 'signalp.1' file to a appropriate location  in your manual
#       system. If you need a compiled version try running:
# 
#            man -d signalp.1 | compress >signalp.Z
# 
#            or:
# 
#            neqn signalp.1 | tbl | nroff -man | col | compress >signalp.Z
# 
#    6. Enjoy ...

SOFTWARE=signalp
VERSION=4.1
ARCHIVE=$SOFTWARE-${VERSION}c.Linux.tar.gz
ARCHIVE_URL=http://www.dropbox.com/s/b4guq6ysyhi7eqm/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd $SOFTWARE_DIR

  # Update signalp program path with dynamic path cmd
  sed -i "s,\$ENV{SIGNALP} = '/usr/cbs/bio/src/$SOFTWARE_DIR',use FindBin;\n    \$ENV{SIGNALP} = \$FindBin::Bin," signalp

  # Update where to store temporary files (writable to all users)
  mkdir tmp
  sed -i "s,my \$outputDir = \"/var/tmp,my \$outputDir = \"\$ENV{SIGNALP}/tmp," signalp

  # Update max number of sequences per run
  sed -i 's/$MAX_ALLOWED_ENTRIES=10000;/$MAX_ALLOWED_ENTRIES=2000000;/' signalp

  # Update Perl script shebang
  sed -i "s,#!/usr/bin/perl,#!/usr/bin/env perl," signalp

  # Test SignalP on the 10 eukaryotic sequences shipped with the package
  set +e
  echo "Testing SignalP..."
  ./signalp -t euk -f short test/euk10.fsa > euk10.fsa.short_out
  ./signalp -t euk -f long test/euk10.fsa > euk10.fsa.long_out
  ./signalp -t euk -f all test/euk10.fsa > euk10.fsa.all_out
  ./signalp -t euk -f summary test/euk10.fsa > euk10.fsa.summary_out
  echo "The output files *_out  should be identical to the corresponding files in '$SOFTWARE_DIR/test', save round number differences"
  echo "diff test/euk10.fsa.short_out euk10.fsa.short_out"
  diff test/euk10.fsa.short_out euk10.fsa.short_out
  echo "diff test/euk10.fsa.long_out euk10.fsa.long_out"
  diff test/euk10.fsa.long_out euk10.fsa.long_out
  echo "diff test/euk10.fsa.all_out euk10.fsa.all_out"
  diff test/euk10.fsa.all_out euk10.fsa.all_out
  echo "diff test/euk10.fsa.summary_out euk10.fsa.summary_out"
  diff test/euk10.fsa.summary_out euk10.fsa.summary_out
  echo "Testing SignalP finished."
  set -e

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
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
