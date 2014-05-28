#!/bin/bash

#
# SignalP
#

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
#  

SOFTWARE=signalp
VERSION=4.1

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE

# Create install directory with permissions if necessary
if [[ ! -d $INSTALL_DIR ]]
then
  mkdir $INSTALL_DIR
  chmod ug+rwX,o+rX $INSTALL_DIR
fi

INSTALL_DOWNLOAD=$INSTALL_DIR/tmp
mkdir $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
ARCHIVE=$SOFTWARE-${VERSION}c.Linux.tar.gz
# If archive was previously downloaded, use the local one, otherwise get it from remote site
if [[ -f ${!INSTALL_HOME}/archive/$ARCHIVE ]]
then
  echo "Archive $ARCHIVE already in ${!INSTALL_HOME}/archive/: using it..."
  cp -a ${!INSTALL_HOME}/archive/$ARCHIVE .
else
  echo "Archive $ARCHIVE not in ${!INSTALL_HOME}/archive/: downloading it..."
  wget http://www.dropbox.com/s/b4guq6ysyhi7eqm/$ARCHIVE -O $ARCHIVE
fi
tar zxvf $ARCHIVE

SOFTWARE_DIR=$SOFTWARE-$VERSION
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

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX,o+rX .
mv -i $SOFTWARE_DIR $INSTALL_DIR
# Store archive if not already present or if different from the previous one
if [[ ! -f ${!INSTALL_HOME}/archive/$ARCHIVE || `diff ${!INSTALL_HOME}/archive/$ARCHIVE $ARCHIVE` ]]
then
  mv -i $ARCHIVE ${!INSTALL_HOME}/archive/
fi

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                \$::env($INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE_DIR
prepend-path    PATH                \$root
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Set module directory path by removing '_INSTALL_HOME' in $INSTALL_HOME and lowercasing the result
MODULE_DIR=${!INSTALL_HOME}/modulefiles/`echo ${INSTALL_HOME/_INSTALL_HOME/} | tr '[:upper:]' '[:lower:]'`/$SOFTWARE

# Create module directory with permissions if necessary
if [[ ! -d $MODULE_DIR ]]
then
  mkdir $MODULE_DIR
  chmod ug+rwX,o+rX $MODULE_DIR
fi

# Add permissions and install module
chmod ug+rwX,o+rX $VERSION .version
mv $VERSION .version $MODULE_DIR

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
