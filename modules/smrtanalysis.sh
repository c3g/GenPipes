#!/bin/bash

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# You should probably also delete this commented-out header and the ## comments
################################################################################


#
# Software_name  memtime.
#

SOFTWARE=smrtanalysis 
VERSION=2.1.1
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
wget http://programs.pacificbiosciences.com/l/1652/2013-11-05/2tqk3c         
./smrtanalysis-2.1.1-centos-5.6.run --extract-only

## I didn't had time to automate this step, but to install smrtanalysis without support for mysql (i.e. smrtportal)
## the following file needs to be modified:
## ./smrtanalysis/install/smrtanalysis-2.1.1.128549/postinstall/bin/conf_main_fab.py
## Comment the following lines:
#------------------------------------
##   #Deal with running services
##      #for svc, pid in util.servicesRunning().iteritems():
##      #    if pid:
##      #        if not promptForKillProcess(svc, pid):
##      #            abort('Unable to end process %s - please run this script with the root account.' % svc)
##
## Then:
##-----------------------------------
##   # createDb()
##
##

## Once done 
export SEYMOUR_HOME=$MUGQIC_INSTALL_HOME/software/$NAME/$NAME-$VERSION/

## Then install by running the installation script. Follow instructions (enter dummy parallel envir. and dummy queues).
## Installation tested with smrtanalysis 2.1.1 and 2.2.0
./smrtanalysis/install/smrtanalysis-2.1.1.128549/etc/scripts/postinstall/configure_smrtanalysis.sh
mv -i ./smrtanalysis/install/smrtanalysis-2.1.1.128549/* $INSTALL_PATH/$SOFTWARE-$VERSION/

# Add permissions and install software
chmod -R 775 $INSTALL_PATH/$SOFTWARE-$VERSION
cd $INSTALL_DOWNLOAD
mv -i $INSTALL_DOWNLOAD/*.run $MUGQIC_INSTALL_HOME/archive

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE-$VERSION \" ; 
}
module-whatis \"$SOFTWARE-$VERSION  \" ;
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION ;
prepend-path    PATH                \$root ;
puts            stderr              \"!!!===> Don't forget to source \${SEYMOUR_HOME}/etc/setup.sh  <===!!!\"  
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

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
