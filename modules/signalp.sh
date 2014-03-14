#!/bin/bash

## NOTE (FL): Tools from cbs.dtu.dk under some sort of danish license and need to be downloaded manually.

SOFTWARE=signalp
VERSION=4.1  ## WARNING: version hard-coded below
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download, extract, build
wget https://www.dropbox.com/s/b4guq6ysyhi7eqm/signalp-4.1c.Linux.tar.gz  -O signalp-4.1c.Linux.tar.gz
tar xvf signalp-4.1c.Linux.tar.gz
cd $SOFTWARE-$VERSION 

# Stupid manual changes required by signalP
sed -i "s,/usr/cbs/bio/src/signalp-$VERSION,$PWD,g" signalp 
mkdir -p tmp
sed -i "s,/var/tmp,$PWD/tmp,g" signalp
sed -i 's/$MAX_ALLOWED_ENTRIES=10000;/$MAX_ALLOWED_ENTRIES=2000000;/g' signalp

./signalp -t euk -f short test/euk10.fsa > euk10.fsa.short_out
./signalp -t euk -f long test/euk10.fsa > euk10.fsa.long_out
./signalp -t euk -f all test/euk10.fsa > euk10.fsa.all_out
 ./signalp -t euk -f summary test/euk10.fsa > euk10.fsa.summary_out



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
# 
cd ..


# Add permissions and install software
chmod -R ug+rwX .
chmod -R o+rX .
mv -i signalp-*.Linux.tar.gz $MUGQIC_INSTALL_HOME/archive 



# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ; 
}
module-whatis \"$SOFTWARE  \" ; 

set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION ;  
prepend-path    PATH                \$root ;
" > $VERSION


################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE


